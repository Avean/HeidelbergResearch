_rd_startup_started_ns = time_ns()
_rd_startup_log = function (label, started_ns)
    elapsed = (time_ns() - started_ns) / 1.0e9
    println("[startup] ", rpad(label, 36, '.'), " ", lpad(string(round(elapsed; digits = 3)), 8), " s")
    flush(stdout)
    return nothing
end

_rd_stage_started_ns = time_ns()
using Pkg
_rd_startup_log("Load Pkg", _rd_stage_started_ns)

_rd_stage_started_ns = time_ns()
Pkg.activate(@__DIR__)
_rd_startup_log("Activate project", _rd_stage_started_ns)

_rd_stage_started_ns = time_ns()
include(joinpath(@__DIR__, "src", "StartupSplash.jl"))
using .StartupSplash
_rd_startup_log("Load startup splash", _rd_stage_started_ns)

_rd_startup_splash = StartupSplash.open_startup_splash!()

# A failure or Cancel in the splash stops the startup here. The splash also
# closes by itself when Julia exits, but not when this file is run again from
# an open REPL, so it is closed explicitly.
try
    # The heavy packages are loaded one at a time before the application
    # modules, which then only bind them, so each gets its own splash stage and
    # timing line. QML goes first: loaded after the plotting and solver packages
    # it invalidates much of their code, and CxxWrap alone took 14 s instead of 2.
    StartupSplash.startup_stage!(_rd_startup_splash, :qml)
    stage_started_ns = time_ns()
    import QML
    _rd_startup_log("Load QML", stage_started_ns)

    StartupSplash.startup_stage!(_rd_startup_splash, :makie)
    stage_started_ns = time_ns()
    import GLMakie, QMLMakie, CairoMakie
    _rd_startup_log("Load Makie", stage_started_ns)

    StartupSplash.startup_stage!(_rd_startup_splash, :solver)
    stage_started_ns = time_ns()
    import DifferentialEquations
    _rd_startup_log("Load DifferentialEquations", stage_started_ns)

    StartupSplash.startup_stage!(_rd_startup_splash, :models)
    stage_started_ns = time_ns()
    include(joinpath(@__DIR__, "src", "ReactionDiffusionApp.jl"))
    _rd_startup_log("Load ReactionDiffusionApp", stage_started_ns)

    StartupSplash.startup_stage!(_rd_startup_splash, :interface)
    stage_started_ns = time_ns()
    include(joinpath(@__DIR__, "src", "QMLInterface.jl"))
    _rd_startup_log("Load ReactionDiffusionQML", stage_started_ns)

    stage_started_ns = time_ns()
    using .ReactionDiffusionQML
    _rd_startup_log("Import QML module", stage_started_ns)

    # The first call compiles run_qml_app and everything it calls.
    StartupSplash.startup_stage!(_rd_startup_splash, :compile)
    Base.invokelatest(
        ReactionDiffusionQML.run_qml_app;
        N = 300,
        startup_started_ns = _rd_startup_started_ns,
        startup_splash = _rd_startup_splash,
        compile_started_ns = time_ns(),
    )
catch err
    err isa StartupSplash.StartupCancelled || rethrow()
    println("[startup] Startup cancelled")
finally
    StartupSplash.close_startup_splash!(_rd_startup_splash)
end
nothing



