module ReactionDiffusionApp

using GLMakie
using DifferentialEquations
using OrdinaryDiffEqSDIRK: TRBDF2
using SparseArrays
using LinearAlgebra
using Random
using Printf
using LaTeXStrings
using Statistics

# ============================================================
# Main module file
# ============================================================

include("Types.jl")
include("Spatial.jl")

# User-friendly model definition layer.
# Must be included before ModelLoader.jl, because model files may use RDModel(...).
include("ModelDSL.jl")

include("ModelLoader.jl")

const MODEL_DIR = normpath(joinpath(@__DIR__, "..", "Models"))
const MODEL_REGISTRY = let started_ns = time_ns()
    registry = load_model_registry(MODEL_DIR)
    elapsed = (time_ns() - started_ns) / 1.0e9
    @printf(
        "[startup] %s %8.3f s\n",
        rpad("Load model registry", 36, '.'),
        elapsed,
    )
    flush(stdout)
    registry
end

include("Simulation.jl")
include("DomainPartition.jl")
include("HeadDetection.jl")
include("HeadConfigurations.jl")
include("SeriesSimulation.jl")

include("PerturbationPanel.jl")
include("PlotPanel.jl")

include("UIRuntime.jl")
include("TopMenu.jl")
include("ControlPanel.jl")
include("PartitionControlPanel.jl")
include("UI.jl")




export run_app

export ModelSpec
export RDModel

export SimulationState
export neumann_laplacian_1d
export periodic_laplacian_1d
export laplacian_1d

export model_files
export create_simulation_state

end
