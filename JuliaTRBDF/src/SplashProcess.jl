# Startup splash process, started by StartupSplash.open_startup_splash!.
#
# Every line on standard input is a stage: start percent, end percent, expected
# seconds and label, separated by tabs. The bar eases from the start towards
# the end of the current stage, and the splash closes when its input ends,
# i.e. when the main process closes the pipe or exits for any reason. Cancel
# exits with the code given as the first argument and, if the second argument
# is a process id rather than 0, first ends that (main) process.

using Observables
using QML

length(ARGS) == 2 || error("Expected the cancel exit code and the process id to end on Cancel.")
const CANCEL_EXIT_CODE = parse(Int, ARGS[1])
const TERMINATE_PID = parse(Int, ARGS[2])
const SPLASH_FILE = normpath(joinpath(@__DIR__, "..", "qml", "Splash.qml"))
isfile(SPLASH_FILE) || error("Startup splash QML file does not exist: $SPLASH_FILE")

# On Windows Qt reads the environment through the C runtime, which does not see
# changes made with ENV; only variables inherited at process start reach it.
# That is why the splash looked different when started from an open REPL.
# qputenv sets the variables where Qt looks.
QML.qputenv("QT_QUICK_CONTROLS_STYLE", QML.QByteArray("Basic"))
QML.qputenv("QSG_RENDER_LOOP", QML.QByteArray("basic"))

mutable struct SplashStage
    start::Float64
    finish::Float64
    seconds::Float64
    received::Float64
end

const stage = SplashStage(0.0, 0.0, 1.0, time())
const stage_text = Observable("Starting…")
const percent = Observable(0.0)
const input_open = Ref(true)
const cancelled = Ref(false)


function read_stage!(line::AbstractString)
    fields = split(line, '\t'; limit = 4)
    length(fields) == 4 || return nothing
    numbers = tryparse.(Float64, fields[1:3])
    any(isnothing, numbers) && return nothing
    stage.start, stage.finish, stage.seconds = numbers
    stage.received = time()
    stage_text[] = fields[4]
    return nothing
end


# Within a stage the bar eases towards its end: about 86 % of the way after the
# expected duration, and never past it, so it keeps moving while the main
# process has nothing new to report.
function eased_percent(now::Float64)
    elapsed = now - stage.received
    progress = min(1 - exp(-2 * elapsed / max(stage.seconds, 0.1)), 0.97)
    return stage.start + (stage.finish - stage.start) * progress
end


function terminate_process(pid::Integer)
    if Sys.iswindows()
        PROCESS_TERMINATE = 0x0001
        handle = ccall(:OpenProcess, stdcall, Ptr{Cvoid}, (UInt32, Cint, UInt32), PROCESS_TERMINATE, 0, pid)
        handle == C_NULL && return nothing
        ccall(:TerminateProcess, stdcall, Cint, (Ptr{Cvoid}, UInt32), handle, 0)
        ccall(:CloseHandle, stdcall, Cint, (Ptr{Cvoid},), handle)
    else
        ccall(:kill, Cint, (Cint, Cint), pid, 9)  # SIGKILL: no crash report
    end
    return nothing
end


reader = @async begin
    try
        for line in eachline(stdin)
            read_stage!(line)
        end
    catch error
        @debug "Splash input failed." exception = error
    end
    input_open[] = false
end

# Stages already waiting in the pipe are read first, and a startup that
# finished while this process was loading never shows the splash at all.
sleep(0.05)
input_open[] || exit(0)

QML.qmlfunction("cancelStartup", () -> (cancelled[] = true; nothing))
QML.loadqml(SPLASH_FILE; startup = QML.JuliaPropertyMap("stage" => stage_text, "percent" => percent))

while input_open[] && !cancelled[]
    QML.process_eventloop_updates()
    QML.process_events()
    value = eased_percent(time())
    value > percent[] && (percent[] = value)
    sleep(0.03)
end

# The main process may be in the middle of a long stage and would only notice
# the Cancel afterwards, so it is ended here when that is allowed.
if cancelled[] && TERMINATE_PID != 0
    println("[startup] Startup cancelled")
    flush(stdout)
    terminate_process(TERMINATE_PID)
end

QML.quit(QML.get_qmlengine())
QML.quit()
QML.cleanup()
QML.process_events()
exit(cancelled[] ? CANCEL_EXIT_CODE : 0)
