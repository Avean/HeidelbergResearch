module ReactionDiffusionApp

using GLMakie
using DifferentialEquations
using SparseArrays
using LinearAlgebra
using Random

# ============================================================
# Main module file
# ============================================================

include("Types.jl")
include("Spatial.jl")
include("ModelLoader.jl")

const MODEL_DIR = normpath(joinpath(@__DIR__, "..", "models"))
const MODEL_REGISTRY = load_model_registry(MODEL_DIR)

include("Simulation.jl")
include("PlotPanel.jl")
include("UI.jl")

export run_app
export ModelSpec
export SimulationState
export neumann_laplacian_1d
export model_files
export create_simulation_state

end