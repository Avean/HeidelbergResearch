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
include(joinpath(@__DIR__, "src", "ReactionDiffusionApp.jl"))
_rd_startup_log("Load ReactionDiffusionApp", _rd_stage_started_ns)

_rd_stage_started_ns = time_ns()
include(joinpath(@__DIR__, "src", "QMLInterface.jl"))
_rd_startup_log("Load ReactionDiffusionQML", _rd_stage_started_ns)

_rd_stage_started_ns = time_ns()
using .ReactionDiffusionQML
_rd_startup_log("Import QML module", _rd_stage_started_ns)

ReactionDiffusionQML.run_qml_app(
    N = 300,
    startup_started_ns = _rd_startup_started_ns,
)
nothing


