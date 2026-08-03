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
    dts = filter(isfinite, [snapshot.dt for snapshot in snapshots])

    return PartitionSnapshot(
        snapshots,
        generation,
        first(snapshots).t,
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
    total_points = total_partition_points(app)
    total_points > 0 || return 0.0

    return partition_base_length(app) *
           app.simulations[segment].N /
           total_points
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
    global_x = reduce(vcat, (sim.x for sim in simulations))

    offsets = cumsum(vcat(0, [sim.N for sim in simulations]))

    for (profile_name, profile_fun) in profiles
        global_values = _evaluate_spatial_profile_for_parameter(
            global_x,
            first_params,
            profile_name,
            profile_fun,
        )

        override_key = spatial_profile_override_key(profile_name)

        for segment in eachindex(simulations)
            indices = (offsets[segment] + 1):offsets[segment + 1]
            simulations[segment].params[override_key] =
                copy(global_values[indices])
        end
    end

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

    left_indices = 1:left_count
    right_indices = (left_count + 1):parent.N
    left_params = slice_partition_params(parent.params, left_indices, parent.N)
    right_params = slice_partition_params(parent.params, right_indices, parent.N)

    left = create_simulation_state_from_data(
        parent.model,
        parent.x[left_indices],
        parent.dx,
        vec(copy(U[left_indices, :])),
        left_params;
        boundary_condition = :neumann,
        displayed_time = displayed_time,
        dtmax = dtmax,
        reltol = reltol,
        abstol = abstol,
    )

    right = create_simulation_state_from_data(
        parent.model,
        parent.x[right_indices],
        parent.dx,
        vec(copy(U[right_indices, :])),
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

    isapprox(left.dx, right.dx; rtol = 1e-10, atol = 0.0) ||
        error("Merged segments must use the same grid spacing.")

    left_U = copy(solution_matrix(left))
    right_U = copy(solution_matrix(right))
    merged_U = vcat(left_U, right_U)
    merged_x = vcat(left.x, right.x)
    displayed_time = current_display_time(left)
    dtmax = min(current_dtmax(left), current_dtmax(right))
    merged_params = merge_partition_params(left, right)

    return create_simulation_state_from_data(
        left.model,
        merged_x,
        left.dx,
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
    app.sim = first(app.simulations)
    refresh_partition_spatial_profile_overrides!(app.simulations)

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
    app.sim = first(app.simulations)
    refresh_partition_spatial_profile_overrides!(app.simulations)

    return nothing
end
