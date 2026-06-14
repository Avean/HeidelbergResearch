# src/Simulation.jl

# ============================================================
# Simulation engine
# ============================================================
#
# This file contains the numerical part of the application.
# It does not create UI elements and does not know how the solution
# is plotted.
#
# Main responsibilities:
#
#     - create the spatial grid,
#     - build the initial condition,
#     - build the ODEProblem,
#     - initialize the time integrator,
#     - step the simulation,
#     - restart the solver after manual perturbations,
#     - control dtmax.
#
# ============================================================


function solution_matrix(sim::SimulationState)
    # Return the current solution reshaped as an N × nvars matrix.
    #
    # Convention:
    #
    #     U[:, 1] = first variable
    #     U[:, 2] = second variable
    #     ...

    return reshape(sim.integrator_ref[].u, sim.N, sim.model.nvars)
end


function solution_matrix(y::AbstractVector, N::Int, nvars::Int)
    # Reshape a raw solver vector as an N × nvars matrix.

    return reshape(y, N, nvars)
end


function make_initial_state(model::ModelSpec, x::Vector{Float64})
    N = length(x)

    params = copy(model.default_params)
    U0 = zeros(Float64, N, model.nvars)

    model.initialize!(U0, x, params)

    return vec(U0), params
end

function make_problem(
    model::ModelSpec,
    y0::Vector{Float64},
    Lap::SparseMatrixCSC{Float64, Int},
    x::Vector{Float64},
    params::Dict{Symbol, Float64},
)
    N = length(x)
    nvars = model.nvars

    function rhs!(dy, y, _p, t)
        U  = reshape(y,  N, nvars)
        dU = reshape(dy, N, nvars)

        model.rhs!(dU, U, Lap, x, params, t)

        return nothing
    end

    return ODEProblem(rhs!, y0, (0.0, Inf))
end

function create_simulation_state(
    model::ModelSpec;
    N::Int = 1000,
    xmin::Float64 = 0.0,
    xmax::Float64 = 1.0,
    dtmax::Float64 = 1e-2,
    reltol::Float64 = 1e-5,
    abstol::Float64 = 1e-7,
)
    # Create a complete simulation state for the selected model.

    validate_model(model)

    x, dx = make_grid_1d(N; xmin = xmin, xmax = xmax)
    Lap = neumann_laplacian_1d(N, dx)

    y0, params = make_initial_state(model, x)

    prob = make_problem(model, y0, Lap, x, params)

    integrator = init(
        prob,
        TRBDF2();
        adaptive = true,
        dt = min(1e-4, dtmax),
        dtmax = dtmax,
        reltol = reltol,
        abstol = abstol,
        save_everystep = false,
    )

    return SimulationState(
        model,
        N,
        x,
        dx,
        Lap,
        params,
        prob,
        Ref{Any}(integrator),
        Ref(0.0),
        Ref(0),
    )
end


function current_internal_dt(sim::SimulationState)
    # Return the current internal/proposed solver time step.
    #
    # If the integrator is temporarily unavailable, return NaN.

    try
        return sim.integrator_ref[].dt
    catch
        return NaN
    end
end


function current_display_time(sim::SimulationState)
    # Return the displayed physical time.
    #
    # This is solver internal time plus accumulated time offset.

    return sim.integrator_ref[].t + sim.time_offset[]
end


function current_dtmax(sim::SimulationState)
    # Return the current maximum allowed solver time step.

    return sim.integrator_ref[].opts.dtmax
end


function set_dtmax!(sim::SimulationState, new_dtmax::Float64)
    # Change the maximum allowed time step of the current integrator.

    isfinite(new_dtmax) && new_dtmax > 0 ||
        error("dtmax must be positive and finite.")

    sim.integrator_ref[].opts.dtmax = new_dtmax

    dt = current_internal_dt(sim)

    if isfinite(dt) && dt > new_dtmax
        set_proposed_dt!(sim.integrator_ref[], new_dtmax)
    end

    return nothing
end


function step_simulation!(sim::SimulationState)
    # Advance the simulation by one accepted solver step.

    step!(sim.integrator_ref[])
    sim.step_counter[] += 1

    shift_time_to_zero_if_needed!(sim)

    return nothing
end


function step_simulation!(sim::SimulationState, nsteps::Int)
    # Advance the simulation by n accepted solver steps.

    nsteps >= 1 || return nothing

    for _ in 1:nsteps
        step_simulation!(sim)
    end

    return nothing
end


function shift_time_to_zero_if_needed!(
    sim::SimulationState;
    threshold::Float64 = 10000.0,
    reltol::Float64 = 1e-5,
    abstol::Float64 = 1e-7,
)
    # Reset internal solver time to zero if it becomes too large.
    #
    # This avoids floating-point precision problems for very large times.
    # The displayed time is preserved by adding the internal solver time
    # to sim.time_offset.

    integrator = sim.integrator_ref[]

    if integrator.t <= threshold
        return nothing
    end

    sim.time_offset[] += integrator.t

    current_u = copy(integrator.u)
    dtmax = current_dtmax(sim)
    dt = min(current_internal_dt(sim), dtmax)

    sim.prob = make_problem(sim.model, current_u, sim.Lap, sim.x, sim.params)

    sim.integrator_ref[] = init(
        sim.prob,
        TRBDF2();
        adaptive = true,
        dt = dt,
        dtmax = dtmax,
        reltol = reltol,
        abstol = abstol,
        save_everystep = false,
    )

    return nothing
end


function restart_after_manual_change!(
    sim::SimulationState,
    ynew::Vector{Float64};
    dt_after_kick::Float64 = 1e-8,
    reltol::Float64 = 1e-5,
    abstol::Float64 = 1e-7,
)
    # Restart the solver after a manual change of the solution.
    #
    # This is important because the adaptive solver should not reuse
    # step-size history from before the perturbation.

    length(ynew) == sim.N * sim.model.nvars ||
        error("New state has wrong length.")

    sim.time_offset[] += sim.integrator_ref[].t

    dtmax = max(dt_after_kick, current_dtmax(sim))

    sim.prob = make_problem(sim.model, copy(ynew), sim.Lap, sim.x, sim.params)

    sim.integrator_ref[] = init(
        sim.prob,
        TRBDF2();
        adaptive = true,
        dt = dt_after_kick,
        dtmax = dtmax,
        reltol = reltol,
        abstol = abstol,
        save_everystep = false,
    )

    sim.step_counter[] = 0

    return nothing
end


function kick_all!(
    sim::SimulationState;
    amount::Float64 = 1.0,
    variable::Int = 1,
)
    # Add a constant perturbation to one selected variable.

    1 <= variable <= sim.model.nvars ||
        error("Invalid variable index.")

    U = copy(solution_matrix(sim))
    U[:, variable] .+= amount

    restart_after_manual_change!(sim, vec(U))

    return nothing
end


function sinusoidal_kick!(
    sim::SimulationState;
    amount::Float64 = 1.0,
    mode::Int = 1,
    variable::Int = 1,
    shift::Float64 = 0.0,
    clip_to_nonnegative::Bool = false,
)
    # Add a sinusoidal perturbation to one selected variable.

    1 <= variable <= sim.model.nvars ||
        error("Invalid variable index.")

    mode >= 1 ||
        error("Mode must be at least 1.")

    U = copy(solution_matrix(sim))

    @. U[:, variable] += amount * sin(2π * mode * sim.x) + shift

    if clip_to_nonnegative
        @. U[:, variable] = max(U[:, variable], 0.0)
    end

    restart_after_manual_change!(sim, vec(U))

    return nothing
end