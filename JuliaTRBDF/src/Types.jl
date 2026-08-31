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

    default_boundary_condition::Union{Nothing, Symbol} = nothing
    # Optional boundary condition selected automatically when the model is
    # chosen. A value of `nothing` preserves the currently selected condition.

    spatial_profile_sets::Vector{Tuple{String, Vector{Tuple{String, Function}}}} =
    Tuple{String, Vector{Tuple{String, Function}}}[]
    # Optional spatial profiles that can be displayed in the UI.

    latex_equations::Vector{String} = String[]
    # Optional LaTeX strings describing the equations, for display in the UI.
end


function validate_model(model::ModelSpec)
    model.nvars >= 1 ||
        error("Model must have at least one variable.")

    length(model.varnames) == model.nvars ||
        error("Length of varnames must be equal to nvars.")

    length(unique(model.varnames)) == length(model.varnames) ||
        error("Variable names must be unique.")

    if !isnothing(model.default_boundary_condition)
        validate_boundary_condition(model.default_boundary_condition)
    end

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

    params::Dict{Symbol, Any}
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


mutable struct PartitionSnapshot
    segments::Vector{SimulationSnapshot}
    generation::Int
    t::Float64
    dt::Float64
    dtmax::Float64
    steps::Int
end


mutable struct SnapshotBuffer
    latest::Base.RefValue{Union{Nothing, PartitionSnapshot}}
    # The newest snapshot produced by the worker thread.

    lock::ReentrantLock
    # Lock protecting access to the latest snapshot.
end


function empty_snapshot_buffer()
    return SnapshotBuffer(
        Ref{Union{Nothing, PartitionSnapshot}}(nothing),
        ReentrantLock(),
    )
end


mutable struct SegmentRuntime
    running::Threads.Atomic{Bool}
    task_ref::Base.RefValue{Union{Nothing, Task}}
    lock::ReentrantLock
    latest_snapshot::Base.RefValue{Union{Nothing, SimulationSnapshot}}
    snapshot_lock::ReentrantLock
end


function empty_segment_runtime()
    return SegmentRuntime(
        Threads.Atomic{Bool}(false),
        Ref{Union{Nothing, Task}}(nothing),
        ReentrantLock(),
        Ref{Union{Nothing, SimulationSnapshot}}(nothing),
        ReentrantLock(),
    )
end


# ============================================================
# In-memory simulation checkpoint
# ============================================================

struct SavedSegmentState
    model_id::Symbol
    x::Vector{Float64}
    dx::Float64
    y::Vector{Float64}
    params::Dict{Symbol, Any}
    boundary_condition::Symbol
    t::Float64
    dt::Float64
    dtmax::Float64
    steps::Int
end


struct SavedSimulationState
    model_id::Symbol
    segments::Vector{SavedSegmentState}
    domain_length_scale::Float64
    requested_dtmax::Float64
    initial_N::Int
    initial_boundary_condition::Symbol
end


# ============================================================
# Plot panel
# ============================================================

mutable struct PlotPanel
    axes::Vector{Axis}
    x_observable::Observable{Vector{Float64}}
    domain_length_scale::Float64
    observables::Vector{Observable{Vector{Float64}}}
    preview_observables::Vector{Observable{Vector{Float64}}}
    perturbation_controls::Vector{Any}
    ui_items::Vector{Any}

    segment_axes::Vector{Vector{Axis}}
    segment_x_observables::Vector{Observable{Vector{Float64}}}
    segment_observables::Vector{Vector{Observable{Vector{Float64}}}}
    segment_preview_observables::Vector{Vector{Observable{Vector{Float64}}}}
    segment_profile_axes::Vector{Vector{Axis}}
    segment_profile_observables::Vector{Vector{Observable{Vector{Float64}}}}
    split_marker_observables::Vector{Observable{Vector{Float64}}}
    split_marker_alpha_observables::Vector{Observable{Float64}}
    split_marker_fade_tokens::Vector{Base.RefValue{Int}}
    segment_status_observables::Vector{Observable{String}}
end

# ============================================================
# Application state
# ============================================================

mutable struct AppState
    sim::SimulationState
    # Current simulation state.

    simulations::Vector{SimulationState}
    # Ordered independent domain segments. app.sim always points to the
    # first segment for shared model and parameter metadata.

    initial_N::Int
    # Number of grid points in the original, unsplit domain.

    initial_boundary_condition::Symbol
    # Boundary condition restored by Reset on the unsplit domain.

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

    segment_runtimes::Vector{SegmentRuntime}
    # One independently scheduled worker runtime per domain segment.

    synchronization_running::Threads.Atomic{Bool}
    synchronization_task_ref::Base.RefValue{Union{Nothing, Task}}
    synchronization_status::Observable{String}

    saved_state::Base.RefValue{Union{Nothing, SavedSimulationState}}
    # One reusable in-memory checkpoint. Model changes invalidate it, while
    # Reset and domain-topology edits intentionally leave it available.

    show_embedded_perturbation_controls::Bool
    # The legacy GLMakie interface keeps its controls inside the plot layout.
    # The QML interface renders the same state in its own bottom drawer.
end
