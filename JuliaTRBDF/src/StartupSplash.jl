module StartupSplash

# Startup splash shown while the application loads. It lives in its own Julia
# process (SplashProcess.jl), so it keeps moving while this process is blocked
# loading packages. Stages go to the splash as lines on its standard input;
# closing that pipe, or this process ending for any reason, closes the splash.
# Cancel comes back as the exit code of the splash process, and when Julia runs
# only this application (a script, not a REPL) the splash also ends this process
# at once, even in the middle of a long stage.

const SPLASH_PROCESS_FILE = normpath(joinpath(@__DIR__, "SplashProcess.jl"))
# Not a code a crash produces (abort() exits with 3 on Windows), so a failing
# splash is never taken for Cancel.
const CANCEL_EXIT_CODE = 42

# Startup stages in order, with their approximate durations in seconds on a
# cold start (averages of full startups on the development machine, with the
# splash process loading at the same time). The bar advances in proportion to
# them, so it moves steadily and never goes back; within a stage the splash
# eases towards the end of the stage.
const STARTUP_STAGES = [
    (key = :qml, label = "Loading Qt and QML", seconds = 8.5),
    (key = :makie, label = "Loading Makie", seconds = 9.5),
    (key = :solver, label = "Loading DifferentialEquations", seconds = 5.5),
    (key = :models, label = "Loading models", seconds = 4.5),
    (key = :interface, label = "Loading application interface", seconds = 1.3),
    (key = :compile, label = "Compiling application code", seconds = 13.0),
    (key = :equations, label = "Rendering model equations", seconds = 16.8),
    (key = :simulation, label = "Creating initial simulation", seconds = 1.6),
    (key = :warm_up, label = "Preparing solver", seconds = 9.2),
    (key = :plots, label = "Preparing plots", seconds = 15.2),
    (key = :connect, label = "Connecting the interface", seconds = 0.3),
    (key = :workspace, label = "Opening workspace", seconds = 1.5),
    (key = :first_frame, label = "Rendering the workspace", seconds = 2.0),
]
const STARTUP_SECONDS = sum(stage.seconds for stage in STARTUP_STAGES)


struct StartupCancelled <: Exception end

Base.showerror(io::IO, ::StartupCancelled) = print(io, "Startup was cancelled in the splash window.")


struct StartupSplashHandle
    process::Base.Process
end


function open_startup_splash!()
    isfile(SPLASH_PROCESS_FILE) ||
        error("Startup splash process file does not exist: $SPLASH_PROCESS_FILE")

    project_root = normpath(joinpath(@__DIR__, ".."))
    # In a REPL, ending the process would close the session, so there Cancel
    # only stops the startup at its next stage (0 = do not terminate).
    terminate_pid = isinteractive() ? 0 : getpid()
    command = `$(Base.julia_cmd()) --project=$project_root --startup-file=no $SPLASH_PROCESS_FILE $CANCEL_EXIT_CODE $terminate_pid`
    # Writes go to the splash's standard input; its output goes to this terminal.
    return StartupSplashHandle(open(command, "w", stdout))
end


startup_cancelled(::Nothing) = false

function startup_cancelled(handle::StartupSplashHandle)
    yield()  # lets the event loop register that the splash process has exited
    return process_exited(handle.process) &&
           handle.process.exitcode == CANCEL_EXIT_CODE
end


startup_stage!(::Nothing, ::Symbol; kwargs...) = nothing

function startup_stage!(
    handle::StartupSplashHandle,
    key::Symbol;
    fraction::Real = 0.0,
    detail::AbstractString = "",
)
    startup_cancelled(handle) && throw(StartupCancelled())

    index = findfirst(stage -> stage.key == key, STARTUP_STAGES)
    index === nothing && error("Unknown startup stage: $key")
    stage = STARTUP_STAGES[index]
    completed = sum((earlier.seconds for earlier in STARTUP_STAGES[1:(index - 1)]); init = 0.0)
    start = completed + clamp(fraction, 0, 1) * stage.seconds
    finish = completed + stage.seconds
    text = stage.label * (isempty(detail) ? "" : " " * detail) * "…"

    line = join(
        (100 * start / STARTUP_SECONDS, 100 * finish / STARTUP_SECONDS, finish - start, text),
        '\t',
    )
    try
        println(handle.process, line)
    catch
        # The splash process has ended (closed or failed); start without it.
    end
    return nothing
end


close_startup_splash!(::Nothing) = nothing

function close_startup_splash!(handle::StartupSplashHandle)
    # The splash closes when its input ends; closing it again is harmless.
    try
        close(handle.process.in)
    catch
    end
    return nothing
end


export StartupCancelled
export StartupSplashHandle
export open_startup_splash!
export startup_stage!
export startup_cancelled
export close_startup_splash!

end
