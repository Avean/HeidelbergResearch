module ReactionDiffusionApp

using GLMakie
using DifferentialEquations
using SparseArrays
using LinearAlgebra
using Random
using Printf
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

const MODEL_DIR = normpath(joinpath(@__DIR__, "..", "models"))
const MODEL_REGISTRY = load_model_registry(MODEL_DIR)

include("Simulation.jl")

include("PerturbationPanel.jl")
include("PlotPanel.jl")

include("UIRuntime.jl")
include("TopMenu.jl")
include("ControlPanel.jl")
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