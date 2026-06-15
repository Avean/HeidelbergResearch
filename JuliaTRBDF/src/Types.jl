# src/Types.jl

# ============================================================
# Model specification
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
    #
    # Signature:
    #
    #     initialize!(U, x, p)
    #
    # where U has size N × nvars.

    rhs!::Function
    # Right-hand side of the semi-discrete ODE system.
    #
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

    boundary_condition::Symbol
    # Boundary condition used by the current simulation.
    # Supported values are currently:
    #
    #     :neumann
    #     :periodic

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
# Snapshot passed from worker thread to UI
# ============================================================

mutable struct SimulationSnapshot
    y::Vector{Float64}
    # Copy of the current solver state.
    # This is safe to read from the UI thread because it is independent
    # of the live integrator.

    N::Int
    # Number of spatial grid points used when the snapshot was created.

    nvars::Int
    # Number of variables in the model used when the snapshot was created.

    model_id::Symbol
    # Identifier of the model used when the snapshot was created.

    generation::Int
    # Application generation number.
    # This prevents old snapshots from a previous model from being drawn
    # after switching models.

    t::Float64
    # Displayed physical simulation time.

    dt::Float64
    # Current internal/proposed solver time step.

    dtmax::Float64
    # Current maximum allowed solver time step.

    steps::Int
    # Number of accepted solver steps since the last full restart.
end


mutable struct SnapshotBuffer
    latest::Base.RefValue{Union{Nothing, SimulationSnapshot}}
    # The newest snapshot produced by the worker thread.

    lock::ReentrantLock
    # Lock protecting access to the latest snapshot.
end


function empty_snapshot_buffer()
    return SnapshotBuffer(
        Ref{Union{Nothing, SimulationSnapshot}}(nothing),
        ReentrantLock(),
    )
end


# ============================================================
# Plot panel
# ============================================================

mutable struct PlotPanel
    axes::Vector{Axis}
    observables::Vector{Observable{Vector{Float64}}}
    preview_observables::Vector{Observable{Vector{Float64}}}
    perturbation_controls::Vector{Any}
    ui_items::Vector{Any}
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
    # Whether the simulation is currently running from the UI perspective.
    # This observable should be updated only from the UI thread.

    dtmax_obs::Observable{Float64}
    # Observable storing the current maximum allowed time step.

    dt_obs::Observable{Float64}
    # Observable storing the currently proposed/internal solver time step.

    time_obs::Observable{Float64}
    # Observable storing the displayed physical simulation time.

    step_counter_obs::Observable{Int}
    # Observable storing the number of solver steps shown in the UI.

    worker_running::Threads.Atomic{Bool}
    # Thread-safe flag used by the simulation worker.

    worker_task_ref::Base.RefValue{Union{Nothing, Task}}
    # Reference to the threaded simulation worker task.

    ui_task_ref::Base.RefValue{Union{Nothing, Task}}
    # Reference to the UI-side snapshot polling task.

    snapshot_buffer::SnapshotBuffer
    # Latest simulation snapshot shared between the worker and UI.

    generation::Threads.Atomic{Int}
    # Incremented whenever the model is switched.
    # Old snapshots with older generations are ignored.

    simlock::ReentrantLock
    # Lock protecting the live solver state.
end