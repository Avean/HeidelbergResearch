# src/Types.jl

# ============================================================
# Model specification
# ============================================================
#
# Every file in models/ must define:
#
#     function build_model()
#         return ModelSpec(...)
#     end
#
# The solver stores the solution as a single vector y.
# Inside model functions we reshape it as:
#
#     U = reshape(y, N, nvars)
#
# Convention:
#
#     U[:, 1] = first variable
#     U[:, 2] = second variable
#     ...
#
# ============================================================

Base.@kwdef struct ModelSpec
    id::Symbol
    # Short internal model identifier, for example :fisher_kpp or :schnakenberg.

    display_name::String
    # Human-readable model name displayed in the user interface.

    nvars::Int
    # Number of dependent variables/equations in the model.

    varnames::Vector{String}
    # Names of dependent variables, for example ["u"] or ["u", "v"].

    default_params::Dict{Symbol, Float64}
    # Default numerical parameters of the model, for example diffusion constants.

    initialize!::Function
    # Function initializing the solution matrix U at t = 0.
    # Signature:
    #
    #     initialize!(U, x, p)
    #
    # where U has size N × nvars.

    rhs!::Function
    # Right-hand side of the semi-discrete ODE system.
    # Signature:
    #
    #     rhs!(dU, U, Lap, x, p, t)
    #
    # where U and dU have size N × nvars.
end


function validate_model(model::ModelSpec)
    model.nvars >= 1 ||
        error("Model must have at least one variable.")

    length(model.varnames) == model.nvars ||
        error("Length of varnames must be equal to nvars.")

    length(unique(model.varnames)) == length(model.varnames) ||
        error("Variable names must be unique.")

    return true
end


# ============================================================
# Simulation state
# ============================================================

mutable struct SimulationState
    model::ModelSpec
    # Currently selected mathematical model.

    N::Int
    # Number of spatial grid points.

    x::Vector{Float64}
    # Spatial grid points.

    dx::Float64
    # Spatial mesh size.

    Lap::SparseMatrixCSC{Float64, Int}
    # Sparse matrix representing the 1D Laplacian.

    params::Dict{Symbol, Float64}
    # Current parameter values used by the model.

    prob::ODEProblem
    # DifferentialEquations.jl ODE problem built from the selected model.

    integrator_ref::Ref{Any}
    # Reference to the current time integrator.
    # A Ref is used because we sometimes replace the integrator after restarting.

    time_offset::Base.RefValue{Float64}
    # Accumulated physical time after internal solver time resets.
    # This avoids loss of floating-point precision for very large times.

    step_counter::Base.RefValue{Int}
    # Number of accepted solver steps since the last full restart.
end


# ============================================================
# Plot panel
# ============================================================

mutable struct PlotPanel
    axes::Vector{Axis}
    # One Makie axis for each model variable.

    observables::Vector{Observable{Vector{Float64}}}
    # Plot data observables, one for each model variable.
    # Updating these observables updates the displayed curves.
end


# ============================================================
# Application state
# ============================================================

mutable struct AppState
    sim::SimulationState
    # Current simulation state.

    plot_panel::PlotPanel
    # Current collection of plots displaying the solution.

    running::Observable{Bool}
    # Whether the simulation is currently running.

    dtmax_obs::Observable{Float64}
    # Observable storing the current maximum allowed time step.

    dt_obs::Observable{Float64}
    # Observable storing the currently proposed/internal solver time step.

    time_obs::Observable{Float64}
    # Observable storing the displayed physical simulation time.

    step_counter_obs::Observable{Int}
    # Observable storing the number of solver steps shown in the UI.

    task_ref::Base.RefValue{Union{Nothing, Task}}
    # Reference to the asynchronous simulation task.
    # It is nothing if the simulation loop has not been started yet.

    simlock::ReentrantLock
    # Lock protecting the solver state from simultaneous UI and simulation updates.
end