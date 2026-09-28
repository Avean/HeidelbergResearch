# src/SeriesSimulation.jl

# ============================================================
# Batch / series simulation core
# ============================================================
#
# This file is deliberately independent from QML.  The QML layer owns the
# controls and threads, while this file owns reproducible perturbations,
# stationary-state checks and one realization of a series.

Base.@kwdef mutable struct SeriesPerturbation
    id::Int
    segment::Int
    variable::Int
    position::Float64
    width_min::Float64
    width_max::Float64
    height_min::Float64
    height_max::Float64
end


Base.@kwdef mutable struct SeriesSettings
    run_count::Int = 100
    # A panel that has not become stationary by maximum_time is reported as
    # not converged.  maximum_steps_per_panel is only a safety net against a
    # stalled solver whose step size collapses before reaching maximum_time.
    maximum_time::Float64 = 1e15
    maximum_steps_per_panel::Int = 10_000_000
    check_interval::Float64 = 1e3
    # Bound on the scaled time derivative max|du/dt| / max(1, max|u|).
    tolerance::Float64 = 1e-6
    required_consecutive_checks::Int = 3
    dtmax::Float64 = 1e10
    seed::UInt64 = 0x3039
    head_variable::Int = 1
    live_preview::Bool = true
    reltol::Float64 = 1e-5
    abstol::Float64 = 1e-7
end


# Temporary diagnostic switch. Enable it from the REPL without reloading the
# application with:
#
#     ReactionDiffusionApp.SERIES_DEBUG_PRINT_RESIDUALS[] = true
#
# Every actual convergence check is printed as one complete line. A dedicated
# lock keeps output from independently running panels from interleaving.
const SERIES_DEBUG_PRINT_RESIDUALS = Ref(false)
const SERIES_DEBUG_PRINT_LOCK = ReentrantLock()


struct SeriesSegmentTemplate
    model::ModelSpec
    x::Vector{Float64}
    dx::Float64
    y::Vector{Float64}
    params::Dict{Symbol, Any}
    boundary_condition::Symbol
    displayed_time::Float64
end


struct SeriesPanelOutcome
    converged::Bool
    cancelled::Bool
    steps::Int
    time::Float64
    residual::Float64
    heads::Vector{DetectedHead}
    final_snapshot::SimulationSnapshot
end


function validate_series_settings(settings::SeriesSettings, nvars::Int)
    settings.run_count >= 1 || error("Number of series runs must be at least one.")
    settings.maximum_time > 0.0 && isfinite(settings.maximum_time) ||
        error("Maximum series time must be positive and finite.")
    settings.maximum_steps_per_panel >= 1 ||
        error("Maximum solver steps per panel must be at least one.")
    settings.check_interval > 0.0 && isfinite(settings.check_interval) ||
        error("Steady-state check interval must be positive and finite.")
    settings.dtmax > 0.0 && isfinite(settings.dtmax) ||
        error("Series dtmax must be positive and finite.")
    settings.tolerance > 0.0 && isfinite(settings.tolerance) ||
        error("Steady-state tolerance must be positive and finite.")
    settings.required_consecutive_checks >= 1 ||
        error("Required consecutive steady-state checks must be at least one.")
    1 <= settings.head_variable <= nvars ||
        error("Head-detection variable is outside the model variable range.")

    return nothing
end


function validate_series_perturbation!(
    perturbation::SeriesPerturbation,
    simulations::AbstractVector{<:SimulationState},
)
    1 <= perturbation.segment <= length(simulations) ||
        error("Perturbation $(perturbation.id) targets an unavailable panel.")
    sim = simulations[perturbation.segment]
    1 <= perturbation.variable <= sim.model.nvars ||
        error("Perturbation $(perturbation.id) targets an unavailable variable.")
    isfinite(perturbation.position) || error("Perturbation position must be finite.")
    isfinite(perturbation.width_min) && isfinite(perturbation.width_max) ||
        error("Perturbation widths must be finite.")
    0.0 < perturbation.width_min <= perturbation.width_max ||
        error("Perturbation width minimum must be positive and not exceed its maximum.")
    perturbation.width_max <= simulation_domain_length(sim) + eps(Float64) ||
        error("Perturbation width exceeds the selected panel length.")
    isfinite(perturbation.height_min) && isfinite(perturbation.height_max) ||
        error("Perturbation heights must be finite.")
    perturbation.height_min <= perturbation.height_max ||
        error("Perturbation height minimum must not exceed its maximum.")

    return nothing
end


function make_series_template(sim::SimulationState)
    return SeriesSegmentTemplate(
        sim.model,
        copy(sim.x),
        sim.dx,
        copy(sim.integrator_ref[].u),
        deepcopy(sim.params),
        sim.boundary_condition,
        current_display_time(sim),
    )
end


function instantiate_series_template(
    template::SeriesSegmentTemplate,
    settings::SeriesSettings,
)
    return create_simulation_state_from_data(
        template.model,
        template.x,
        template.dx,
        template.y,
        deepcopy(template.params);
        boundary_condition = template.boundary_condition,
        displayed_time = template.displayed_time,
        dtmax = settings.dtmax,
        reltol = settings.reltol,
        abstol = settings.abstol,
    )
end


function _uniform_between(rng::AbstractRNG, left::Float64, right::Float64)
    return left + (right - left) * rand(rng)
end


function apply_series_perturbations!(
    simulations::AbstractVector{<:SimulationState},
    perturbations::AbstractVector{<:SeriesPerturbation},
    rng::AbstractRNG,
    settings::SeriesSettings,
)
    isempty(perturbations) && error("A series requires at least one perturbation.")
    ynew = [copy(sim.integrator_ref[].u) for sim in simulations]

    for perturbation in perturbations
        validate_series_perturbation!(perturbation, simulations)
        sim = simulations[perturbation.segment]
        width = _uniform_between(rng, perturbation.width_min, perturbation.width_max)
        height = _uniform_between(rng, perturbation.height_min, perturbation.height_max)
        mask = local_perturbation_mask(
            sim.x,
            perturbation.position,
            width;
            boundary_condition = sim.boundary_condition,
        )
        U = reshape(ynew[perturbation.segment], sim.N, sim.model.nvars)

        for index in eachindex(mask)
            mask[index] || continue
            U[index, perturbation.variable] += height * rand(rng)
        end
    end

    for segment in eachindex(simulations)
        restart_after_manual_change!(
            simulations[segment],
            ynew[segment];
            reltol = settings.reltol,
            abstol = settings.abstol,
        )
        simulations[segment].step_counter[] = 0
    end

    return nothing
end


function series_stationarity_residual(sim::SimulationState)
    integrator = sim.integrator_ref[]
    du = similar(integrator.u)
    integrator.f(du, integrator.u, integrator.p, integrator.t)
    U = reshape(integrator.u, sim.N, sim.model.nvars)
    dU = reshape(du, sim.N, sim.model.nvars)
    maximum_rate = 0.0

    # The time derivative vanishes exactly at a stationary state, so this test
    # does not depend on the adaptive step size or on the check interval, and
    # a periodic solution cannot pass it by being sampled at its period.
    # Normalising each variable separately ensures that a high-amplitude
    # variable cannot hide continued evolution in a small-amplitude one.
    for variable in 1:sim.model.nvars
        scale = max(1.0, maximum(abs, @view(U[:, variable])))
        maximum_rate = max(maximum_rate, maximum(abs, @view(dU[:, variable])) / scale)
    end

    return maximum_rate
end


function step_series_panel!(sim::SimulationState)
    # Unlike step_simulation!, never shift the solver time back to zero: the
    # panel loop relies on integrator.t for check intervals and maximum_time.
    step!(sim.integrator_ref[])
    sim.step_counter[] += 1

    return nothing
end


function debug_series_residual(
    run_index::Int,
    segment::Int,
    kind::Symbol,
    sim::SimulationState,
    residual::Float64,
    stable_checks::Int,
    settings::SeriesSettings,
)
    SERIES_DEBUG_PRINT_RESIDUALS[] || return nothing

    lock(SERIES_DEBUG_PRINT_LOCK)
    try
        @printf(
            "[series R] run=%d panel=%d kind=%s t=%.16g step=%d R=%.16e tolerance=%.16e passed=%d/%d\n",
            run_index,
            segment,
            String(kind),
            sim.integrator_ref[].t,
            sim.step_counter[],
            residual,
            settings.tolerance,
            stable_checks,
            settings.required_consecutive_checks,
        )
        flush(stdout)
    finally
        unlock(SERIES_DEBUG_PRINT_LOCK)
    end

    return nothing
end


function run_series_panel!(
    sim::SimulationState,
    settings::SeriesSettings,
    generation::Int;
    run_index::Int = 0,
    segment::Int = 0,
    cancelled::Function = () -> false,
    on_snapshot::Function = _ -> nothing,
    on_residual::Function = _ -> nothing,
)
    validate_series_settings(settings, sim.model.nvars)
    stable_checks = 0
    integrator = sim.integrator_ref[]
    converged = false
    residual = NaN
    last_check_time = integrator.t
    last_residual_time = NaN
    last_snapshot_report_ns = UInt64(0)

    while integrator.t < settings.maximum_time
        cancelled() && break
        sim.step_counter[] >= settings.maximum_steps_per_panel && break
        step_series_panel!(sim)
        settings.live_preview && on_snapshot(make_snapshot(sim, generation))

        now_ns = time_ns()
        if !settings.live_preview &&
           now_ns - last_snapshot_report_ns >= UInt64(500_000_000)
            on_snapshot(make_snapshot(sim, generation))
            last_snapshot_report_ns = now_ns
        end

        # check_interval is a minimum separation between convergence checks,
        # not a solver stop. A single accepted step may cross many intervals;
        # it still produces exactly one check at the accepted endpoint.
        if integrator.t - last_check_time >= settings.check_interval
            residual = series_stationarity_residual(sim)
            stable_checks = residual < settings.tolerance ? stable_checks + 1 : 0
            last_check_time = integrator.t
            last_residual_time = integrator.t
            settings.live_preview || on_snapshot(make_snapshot(sim, generation))
            on_residual(residual)
            debug_series_residual(
                run_index,
                segment,
                :check,
                sim,
                residual,
                stable_checks,
                settings,
            )

            if stable_checks >= settings.required_consecutive_checks
                converged = true
                break
            end
        end
    end

    # If the final accepted state was not already checked, report its residual
    # for diagnostics without counting it as another convergence check.
    if last_residual_time != integrator.t
        residual = series_stationarity_residual(sim)
        on_residual(residual)
        debug_series_residual(
            run_index,
            segment,
            :final,
            sim,
            residual,
            stable_checks,
            settings,
        )
    end
    heads = DetectedHead[]

    if converged
        U = solution_matrix(sim)
        heads = detect_heads(
            @view(U[:, settings.head_variable]),
            sim.x;
            boundary_condition = sim.boundary_condition,
        )
    end

    snapshot = make_snapshot(sim, generation)
    on_snapshot(snapshot)
    return SeriesPanelOutcome(
        converged,
        !converged && cancelled(),
        sim.step_counter[],
        integrator.t,
        residual,
        heads,
        snapshot,
    )
end


function run_series_realization!(
    templates::AbstractVector{<:SeriesSegmentTemplate},
    perturbations::AbstractVector{<:SeriesPerturbation},
    settings::SeriesSettings,
    generation::Int,
    run_index::Int;
    cancelled::Function = () -> false,
    on_snapshot::Function = (_, _) -> nothing,
    on_residual::Function = (_, _) -> nothing,
)
    isempty(templates) && error("A series realization needs at least one panel.")
    validate_series_settings(settings, first(templates).model.nvars)
    run_rng = Xoshiro(settings.seed + UInt64(run_index - 1))
    simulations = [instantiate_series_template(template, settings) for template in templates]
    apply_series_perturbations!(simulations, perturbations, run_rng, settings)
    tasks = Task[]

    # Publish a consistent perturbed initial condition before any panel starts
    # advancing.  This keeps the optional live preview from mixing the last
    # realization with the next one.
    for segment in eachindex(simulations)
        on_snapshot(segment, make_snapshot(simulations[segment], generation))
    end

    for segment in eachindex(simulations)
        task = Threads.@spawn run_series_panel!(
            simulations[segment],
            settings,
            generation;
            run_index = run_index,
            segment = segment,
            cancelled = cancelled,
            on_snapshot = snapshot -> on_snapshot(segment, snapshot),
            on_residual = residual -> on_residual(segment, residual),
        )
        push!(tasks, task)
    end

    outcomes = SeriesPanelOutcome[]
    for task in tasks
        push!(outcomes, fetch(task))
    end

    return outcomes
end
