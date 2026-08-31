# src/DomainPartition.jl

# ============================================================
# Partitioned one-dimensional domains
# ============================================================


function total_partition_points(app::AppState)
    return sum(sim.N for sim in app.simulations)
end


function make_partition_snapshot(
    simulations::Vector{SimulationState},
    generation::Int,
)
    snapshots = [make_snapshot(sim, generation) for sim in simulations]
    return partition_snapshot_from_segments(snapshots, generation)
end


function partition_snapshot_from_segments(
    snapshots::Vector{SimulationSnapshot},
    generation::Int,
)
    isempty(snapshots) && error("A partition snapshot needs at least one segment.")
    dts = filter(isfinite, [snapshot.dt for snapshot in snapshots])

    return PartitionSnapshot(
        snapshots,
        generation,
        maximum(snapshot.t for snapshot in snapshots),
        isempty(dts) ? NaN : minimum(dts),
        minimum(snapshot.dtmax for snapshot in snapshots),
        maximum(snapshot.steps for snapshot in snapshots),
    )
end


function step_partition_synchronized!(
    simulations::Vector{SimulationState},
    nsteps::Int,
)
    isempty(simulations) && return nothing
    nsteps >= 1 || return nothing

    if length(simulations) == 1
        step_simulation!(first(simulations), nsteps)
        return nothing
    end

    for _ in 1:nsteps
        proposed_dts = Float64[]

        for sim in simulations
            dt = abs(current_internal_dt(sim))

            if isfinite(dt) && dt > 0
                push!(proposed_dts, min(dt, current_dtmax(sim)))
            end
        end

        isempty(proposed_dts) &&
            error("No valid time step is available for the partition.")

        common_dt = minimum(proposed_dts)

        for sim in simulations
            step!(sim.integrator_ref[], common_dt, true)
            sim.step_counter[] += 1
            shift_time_to_zero_if_needed!(sim)
        end
    end

    return nothing
end


function partition_base_length(app::AppState)
    first_sim = first(app.simulations)
    last_sim = last(app.simulations)
    xmin = first(first_sim.x)
    xmax = last(last_sim.x)

    if length(app.simulations) == 1 &&
       first_sim.boundary_condition == :periodic
        return (xmax - xmin) + first_sim.dx
    end

    return xmax - xmin
end


function segment_base_length(app::AppState, segment::Int)
    sim = app.simulations[segment]

    if length(app.simulations) == 1 && sim.boundary_condition == :periodic
        return (last(sim.x) - first(sim.x)) + sim.dx
    end

    return last(sim.x) - first(sim.x)
end


function segment_physical_length(sim::SimulationState)
    interval_count = sim.boundary_condition == :periodic ? sim.N : sim.N - 1
    return interval_count * sim.dx
end


function rebase_partition_coordinates!(
    simulations::Vector{SimulationState};
    origin::Float64 = 0.0,
)
    cursor = origin

    for sim in simulations
        sim.N == length(sim.x) ||
            error("Segment grid size does not match its coordinate vector.")
        isfinite(sim.dx) && sim.dx > 0 ||
            error("Segment grid spacing must be positive and finite.")

        # Mutate the existing vector rather than replacing it: the ODE RHS
        # closure holds this same x object and immediately sees the new origin.
        for index in eachindex(sim.x)
            sim.x[index] = cursor + (index - 1) * sim.dx
        end

        cursor += segment_physical_length(sim)
    end

    return nothing
end


function segment_display_coordinates(
    app::AppState,
    segment::Int,
    domain_length_scale::Real,
)
    sim = app.simulations[segment]
    displayed_length =
        segment_base_length(app, segment) * Float64(domain_length_scale)

    if sim.boundary_condition == :periodic && length(app.simulations) == 1
        step = displayed_length / sim.N
        return collect(range(0.0; step = step, length = sim.N))
    end

    return collect(range(0.0, displayed_length; length = sim.N))
end


function _is_spatial_profile_override_key(key::Symbol)
    return startswith(String(key), SPATIAL_PROFILE_OVERRIDE_PREFIX)
end


function clear_spatial_profile_overrides!(sim::SimulationState)
    for key in collect(keys(sim.params))
        if _is_spatial_profile_override_key(key)
            delete!(sim.params, key)
        end
    end

    return nothing
end


function slice_partition_params(
    params::AbstractDict{Symbol},
    indices,
    parent_N::Int,
)
    child_params = Dict{Symbol, Any}(params)

    for (key, value) in params
        _is_spatial_profile_override_key(key) || continue
        profile_values = Float64.(collect(value))
        length(profile_values) == parent_N ||
            error("Spatial profile override has wrong length before split.")
        child_params[key] = copy(profile_values[indices])
    end

    return child_params
end


function merge_partition_params(
    left::SimulationState,
    right::SimulationState,
)
    merged_params = Dict{Symbol, Any}(left.params)
    override_keys = union(
        filter(_is_spatial_profile_override_key, keys(left.params)),
        filter(_is_spatial_profile_override_key, keys(right.params)),
    )

    for key in override_keys
        haskey(left.params, key) && haskey(right.params, key) ||
            error("Both merged segments must contain the same spatial profile overrides.")

        left_values = Float64.(collect(left.params[key]))
        right_values = Float64.(collect(right.params[key]))
        length(left_values) == left.N ||
            error("Left spatial profile override has wrong length before merge.")
        length(right_values) == right.N ||
            error("Right spatial profile override has wrong length before merge.")

        merged_params[key] = vcat(left_values, right_values)
    end

    return merged_params
end


function interpolate_partition_values(
    source_x::AbstractVector{<:Real},
    source_values::AbstractVector{<:Real},
    target_x::AbstractVector{<:Real};
    periodic::Bool = false,
    period::Union{Nothing, Float64} = nothing,
)
    length(source_x) == length(source_values) ||
        error("Interpolation source coordinates and values must have equal lengths.")
    length(source_x) >= 2 ||
        error("At least two source points are required for interpolation.")

    x_values = Float64.(collect(source_x))
    y_values = Float64.(collect(source_values))

    if periodic
        period_value = isnothing(period) ?
            (last(x_values) - first(x_values)) + (x_values[2] - x_values[1]) :
            period
        x_values = vcat(x_values, first(x_values) + period_value)
        y_values = vcat(y_values, first(y_values))
    end

    result = Vector{Float64}(undef, length(target_x))

    for (target_index, raw_x) in enumerate(target_x)
        x = Float64(raw_x)

        if x <= first(x_values)
            result[target_index] = first(y_values)
            continue
        elseif x >= last(x_values)
            result[target_index] = last(y_values)
            continue
        end

        right_index = searchsortedfirst(x_values, x)
        left_index = right_index - 1
        x0 = x_values[left_index]
        x1 = x_values[right_index]
        weight = (x - x0) / (x1 - x0)
        result[target_index] =
            (1.0 - weight) * y_values[left_index] +
            weight * y_values[right_index]
    end

    return result
end


function resample_partition_params(
    params::AbstractDict{Symbol},
    source_x::AbstractVector{<:Real},
    target_x::AbstractVector{<:Real};
    periodic::Bool = false,
    period::Union{Nothing, Float64} = nothing,
)
    target_params = Dict{Symbol, Any}(params)

    for (key, value) in params
        _is_spatial_profile_override_key(key) || continue
        target_params[key] = interpolate_partition_values(
            source_x,
            Float64.(collect(value)),
            target_x;
            periodic = periodic,
            period = period,
        )
    end

    return target_params
end


function refresh_partition_spatial_profile_overrides!(
    simulations::Vector{SimulationState},
)
    isempty(simulations) && return nothing

    model = first(simulations).model

    for sim in simulations
        sim.model.id == model.id ||
            error("All domain segments must use the same model.")

        clear_spatial_profile_overrides!(sim)
    end

    isempty(model.spatial_profile_sets) &&
        return nothing

    first_params = first(simulations).params
    set_index = _active_spatial_profile_set_index(
        first_params,
        model.spatial_profile_sets,
    )
    _, profiles = model.spatial_profile_sets[set_index]
    reference_N = first(simulations).N
    global_xmin = first(first(simulations).x)
    global_xmax = last(last(simulations).x)
    global_x = collect(range(global_xmin, global_xmax; length = reference_N))

    for (profile_name, profile_fun) in profiles
        global_values = _evaluate_spatial_profile_for_parameter(
            global_x,
            first_params,
            profile_name,
            profile_fun,
        )

        override_key = spatial_profile_override_key(profile_name)

        for sim in simulations
            sim.params[override_key] = interpolate_partition_values(
                global_x,
                global_values,
                sim.x,
            )
        end
    end

    return nothing
end


function partition_spatial_profile_overrides_are_valid(
    simulations::Vector{SimulationState},
)
    isempty(simulations) && return true
    model = first(simulations).model
    isempty(model.spatial_profile_sets) && return true

    first_params = first(simulations).params
    set_index = _active_spatial_profile_set_index(
        first_params,
        model.spatial_profile_sets,
    )
    _, profiles = model.spatial_profile_sets[set_index]

    for sim in simulations
        sim.model.id == model.id || return false
        _active_spatial_profile_set_index(
            sim.params,
            sim.model.spatial_profile_sets,
        ) == set_index || return false

        for (profile_name, _) in profiles
            override_key = spatial_profile_override_key(profile_name)
            haskey(sim.params, override_key) || return false

            values = sim.params[override_key]
            values isa AbstractVector || return false
            length(values) == sim.N || return false
        end
    end

    return true
end


function ensure_partition_spatial_profile_overrides!(
    simulations::Vector{SimulationState},
)
    partition_spatial_profile_overrides_are_valid(simulations) ||
        refresh_partition_spatial_profile_overrides!(simulations)

    return nothing
end


function split_simulation_state(
    parent::SimulationState,
    left_count::Int;
    reltol::Float64 = 1e-5,
    abstol::Float64 = 1e-7,
)
    2 <= left_count <= parent.N - 2 ||
        error("Both child domains must contain at least two points.")

    U = copy(solution_matrix(parent))
    displayed_time = current_display_time(parent)
    dtmax = current_dtmax(parent)
    parent_periodic = parent.boundary_condition == :periodic
    parent_length = parent_periodic ?
        (last(parent.x) - first(parent.x)) + parent.dx :
        last(parent.x) - first(parent.x)
    split_fraction = left_count / parent.N
    split_x = first(parent.x) + split_fraction * parent_length
    parent_right = first(parent.x) + parent_length

    left_x = collect(range(first(parent.x), split_x; length = parent.N))
    right_x = collect(range(split_x, parent_right; length = parent.N))
    left_U = zeros(Float64, parent.N, parent.model.nvars)
    right_U = similar(left_U)

    for variable in 1:parent.model.nvars
        left_U[:, variable] = interpolate_partition_values(
            parent.x,
            U[:, variable],
            left_x;
            periodic = parent_periodic,
            period = parent_length,
        )
        right_U[:, variable] = interpolate_partition_values(
            parent.x,
            U[:, variable],
            right_x;
            periodic = parent_periodic,
            period = parent_length,
        )
    end

    left_params = resample_partition_params(
        parent.params,
        parent.x,
        left_x;
        periodic = parent_periodic,
        period = parent_length,
    )
    right_params = resample_partition_params(
        parent.params,
        parent.x,
        right_x;
        periodic = parent_periodic,
        period = parent_length,
    )
    left_dx = (last(left_x) - first(left_x)) / (parent.N - 1)
    right_dx = (last(right_x) - first(right_x)) / (parent.N - 1)

    left = create_simulation_state_from_data(
        parent.model,
        left_x,
        left_dx,
        vec(left_U),
        left_params;
        boundary_condition = :neumann,
        displayed_time = displayed_time,
        dtmax = dtmax,
        reltol = reltol,
        abstol = abstol,
    )

    right = create_simulation_state_from_data(
        parent.model,
        right_x,
        right_dx,
        vec(right_U),
        right_params;
        boundary_condition = :neumann,
        displayed_time = displayed_time,
        dtmax = dtmax,
        reltol = reltol,
        abstol = abstol,
    )

    return left, right
end


function merge_simulation_states(
    left::SimulationState,
    right::SimulationState;
    boundary_condition::Symbol = :neumann,
    reltol::Float64 = 1e-5,
    abstol::Float64 = 1e-7,
)
    left.model.id == right.model.id ||
        error("Only segments using the same model can be merged.")

    left_U = copy(solution_matrix(left))
    right_U = copy(solution_matrix(right))
    target_N = max(left.N, right.N)
    merged_start = first(left.x)
    left_length = segment_physical_length(left)
    merged_length = left_length + segment_physical_length(right)
    interface_x = merged_start + left_length

    if boundary_condition == :periodic
        merged_dx = merged_length / target_N
        merged_x = collect(
            range(merged_start; step = merged_dx, length = target_N),
        )
    else
        merged_dx = merged_length / (target_N - 1)
        merged_x = collect(
            range(merged_start; step = merged_dx, length = target_N),
        )
    end
    merged_U = zeros(Float64, target_N, left.model.nvars)

    for variable in 1:left.model.nvars
        for index in eachindex(merged_x)
            if merged_x[index] <= interface_x
                merged_U[index, variable] = interpolate_partition_values(
                    left.x,
                    left_U[:, variable],
                    [merged_x[index]],
                )[1]
            else
                merged_U[index, variable] = interpolate_partition_values(
                    right.x,
                    right_U[:, variable],
                    [merged_x[index]],
                )[1]
            end
        end
    end

    displayed_time = current_display_time(left)
    dtmax = min(current_dtmax(left), current_dtmax(right))
    merged_params = Dict{Symbol, Any}(left.params)

    for key in union(
        filter(_is_spatial_profile_override_key, keys(left.params)),
        filter(_is_spatial_profile_override_key, keys(right.params)),
    )
        left_values = Float64.(collect(left.params[key]))
        right_values = Float64.(collect(right.params[key]))
        values = Vector{Float64}(undef, target_N)

        for index in eachindex(merged_x)
            if merged_x[index] <= interface_x
                values[index] = interpolate_partition_values(
                    left.x,
                    left_values,
                    [merged_x[index]],
                )[1]
            else
                values[index] = interpolate_partition_values(
                    right.x,
                    right_values,
                    [merged_x[index]],
                )[1]
            end
        end

        merged_params[key] = values
    end
    return create_simulation_state_from_data(
        left.model,
        merged_x,
        merged_dx,
        vec(merged_U),
        merged_params;
        boundary_condition = boundary_condition,
        displayed_time = displayed_time,
        dtmax = dtmax,
        reltol = reltol,
        abstol = abstol,
    )
end


function split_domain_segment!(
    app::AppState,
    segment::Int,
    left_count::Int,
)
    1 <= segment <= length(app.simulations) ||
        error("Invalid segment index.")

    left, right = split_simulation_state(
        app.simulations[segment],
        left_count,
    )

    splice!(app.simulations, segment:segment, (left, right))
    splice!(
        app.segment_runtimes,
        segment:segment,
        (empty_segment_runtime(), empty_segment_runtime()),
    )
    app.sim = first(app.simulations)

    return nothing
end


function merge_domain_segments!(
    app::AppState,
    left_segment::Int,
)
    1 <= left_segment < length(app.simulations) ||
        error("Invalid merge boundary.")

    merged_boundary_condition =
        length(app.simulations) == 2 ?
        app.initial_boundary_condition :
        :neumann

    merged = merge_simulation_states(
        app.simulations[left_segment],
        app.simulations[left_segment + 1];
        boundary_condition = merged_boundary_condition,
    )

    splice!(
        app.simulations,
        left_segment:(left_segment + 1),
        (merged,),
    )
    splice!(
        app.segment_runtimes,
        left_segment:(left_segment + 1),
        (empty_segment_runtime(),),
    )
    app.sim = first(app.simulations)

    return nothing
end


function swap_adjacent_domain_segments!(
    app::AppState,
    left_segment::Int,
)
    1 <= left_segment < length(app.simulations) ||
        error("Invalid swap boundary.")
    length(app.simulations) == length(app.segment_runtimes) ||
        error("Segment runtime count does not match the simulation partition.")

    right_segment = left_segment + 1
    app.simulations[left_segment], app.simulations[right_segment] =
        app.simulations[right_segment], app.simulations[left_segment]
    app.segment_runtimes[left_segment], app.segment_runtimes[right_segment] =
        app.segment_runtimes[right_segment], app.segment_runtimes[left_segment]

    rebase_partition_coordinates!(app.simulations)
    app.sim = first(app.simulations)

    return nothing
end


function delete_domain_segment!(app::AppState, segment::Int)
    length(app.simulations) > 1 ||
        error("The last remaining domain segment cannot be deleted.")
    1 <= segment <= length(app.simulations) ||
        error("Invalid segment index.")
    length(app.simulations) == length(app.segment_runtimes) ||
        error("Segment runtime count does not match the simulation partition.")

    deleteat!(app.simulations, segment)
    deleteat!(app.segment_runtimes, segment)
    rebase_partition_coordinates!(app.simulations)
    app.sim = first(app.simulations)

    return nothing
end
