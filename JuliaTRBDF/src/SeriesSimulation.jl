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
    maximum_steps_per_panel::Int = 1000
    check_interval::Float64 = 1.0
    tolerance::Float64 = 1e-8
    required_consecutive_checks::Int = 3
    seed::UInt64 = 0x3039
    head_variable::Int = 1
    live_preview::Bool = false
    preview_interval_seconds::Float64 = 0.5
    reltol::Float64 = 1e-5
    abstol::Float64 = 1e-7
end


struct SeriesSegmentTemplate
    model::ModelSpec
    x::Vector{Float64}
    dx::Float64
    y::Vector{Float64}
    params::Dict{Symbol, Any}
    boundary_condition::Symbol
    displayed_time::Float64
    dtmax::Float64
end


struct SeriesPanelOutcome
    converged::Bool
    cancelled::Bool
    steps::Int
    heads::Vector{DetectedHead}
    final_snapshot::SimulationSnapshot
end


function validate_series_settings(settings::SeriesSettings, nvars::Int)
    settings.run_count >= 1 || error("Number of series runs must be at least one.")
    settings.maximum_steps_per_panel >= 1 ||
        error("Maximum solver steps per panel must be at least one.")
    settings.check_interval > 0.0 && isfinite(settings.check_interval) ||
        error("Steady-state check interval must be positive and finite.")
    settings.tolerance > 0.0 && isfinite(settings.tolerance) ||
        error("Steady-state tolerance must be positive and finite.")
    settings.required_consecutive_checks >= 1 ||
        error("Required consecutive steady-state checks must be at least one.")
    settings.preview_interval_seconds > 0.0 ||
        error("Live-preview interval must be positive.")
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
    dtmax = current_dtmax(sim)
    isfinite(dtmax) && dtmax > 0.0 || error("Cannot capture a series template with invalid dtmax.")

    return SeriesSegmentTemplate(
        sim.model,
        copy(sim.x),
        sim.dx,
        copy(sim.integrator_ref[].u),
        deepcopy(sim.params),
        sim.boundary_condition,
        current_display_time(sim),
        dtmax,
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
        dtmax = template.dtmax,
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
        restart_after_manual_change!(simulations[segment], ynew[segment])
        simulations[segment].step_counter[] = 0
    end

    return nothing
end


function series_relative_change(
    previous::AbstractVector{<:Real},
    current::AbstractVector{<:Real},
    N::Int,
    nvars::Int,
)
    length(previous) == length(current) == N * nvars ||
        error("Series steady-state vectors have incompatible lengths.")
    previous_U = reshape(previous, N, nvars)
    current_U = reshape(current, N, nvars)
    maximum_change = 0.0

    # Normalising each variable separately ensures that a high-amplitude
    # variable cannot hide continued evolution in a small-amplitude one.
    for variable in 1:nvars
        previous_column = @view previous_U[:, variable]
        current_column = @view current_U[:, variable]
        scale = max(
            1.0,
            maximum(abs, previous_column),
            maximum(abs, current_column),
        )
        difference = current_column .- previous_column
        maximum_change = max(
            maximum_change,
            maximum(abs, difference) / scale,
        )
    end

    return maximum_change
end


function run_series_panel!(
    sim::SimulationState,
    settings::SeriesSettings,
    generation::Int;
    cancelled::Function = () -> false,
    on_snapshot::Function = _ -> nothing,
)
    validate_series_settings(settings, sim.model.nvars)
    previous = copy(sim.integrator_ref[].u)
    stable_checks = 0
    last_preview_ns = time_ns() - round(UInt64, settings.preview_interval_seconds * 1e9)
    integrator = sim.integrator_ref[]

    while sim.step_counter[] < settings.maximum_steps_per_panel
        cancelled() && break
        checkpoint = integrator.t + settings.check_interval
        add_tstop!(integrator, checkpoint)

        while integrator.t < checkpoint
            cancelled() && break
            sim.step_counter[] >= settings.maximum_steps_per_panel && break
            step_simulation!(sim)

            if settings.live_preview
                now_ns = time_ns()
                if now_ns - last_preview_ns >= round(UInt64, settings.preview_interval_seconds * 1e9)
                    on_snapshot(make_snapshot(sim, generation))
                    last_preview_ns = now_ns
                end
            end
        end

        cancelled() && break
        sim.step_counter[] >= settings.maximum_steps_per_panel && break

        current = copy(integrator.u)
        change = series_relative_change(previous, current, sim.N, sim.model.nvars)
        stable_checks = change < settings.tolerance ? stable_checks + 1 : 0
        previous = current

        if stable_checks >= settings.required_consecutive_checks
            U = solution_matrix(sim)
            heads = detect_heads(
                @view(U[:, settings.head_variable]),
                sim.x;
                boundary_condition = sim.boundary_condition,
            )
            snapshot = make_snapshot(sim, generation)
            on_snapshot(snapshot)
            return SeriesPanelOutcome(true, false, sim.step_counter[], heads, snapshot)
        end
    end

    snapshot = make_snapshot(sim, generation)
    on_snapshot(snapshot)
    return SeriesPanelOutcome(
        false,
        cancelled(),
        sim.step_counter[],
        DetectedHead[],
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
)
    isempty(templates) && error("A series realization needs at least one panel.")
    validate_series_settings(settings, first(templates).model.nvars)
    run_rng = Xoshiro(settings.seed + UInt64(run_index - 1))
    simulations = [instantiate_series_template(template, settings) for template in templates]
    apply_series_perturbations!(simulations, perturbations, run_rng)
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
            cancelled = cancelled,
            on_snapshot = snapshot -> on_snapshot(segment, snapshot),
        )
        push!(tasks, task)
    end

    outcomes = SeriesPanelOutcome[]
    for task in tasks
        push!(outcomes, fetch(task))
    end

    return outcomes
end
