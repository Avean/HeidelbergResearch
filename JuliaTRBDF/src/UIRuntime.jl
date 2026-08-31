# src/UIRuntime.jl

# ============================================================
# UI runtime logic
# ============================================================
#
# This file contains non-layout UI logic:
#
#     - snapshot-based refresh,
#     - worker thread control,
#     - safe simulation mutations,
#     - model switching,
#     - boundary-condition switching,
#     - reset and constant initial conditions,
#     - local perturbation application,
#     - active spatial profile switching.
#
# It should not create buttons, menus, or axes.
#
# ============================================================


# ============================================================
# Snapshot-based UI refresh
# ============================================================

function snapshot_matches_current_app(
    app::AppState,
    snapshot::PartitionSnapshot,
)
    snapshot.generation == app.generation[] || return false
    length(snapshot.segments) == length(app.simulations) || return false
    length(app.plot_panel.segment_observables) == length(app.simulations) ||
        return false
    length(app.plot_panel.segment_status_observables) == length(app.simulations) ||
        return false

    return all(
        segment_snapshot.generation == app.generation[] &&
        segment_snapshot.model_id == app.sim.model.id &&
        segment_snapshot.N == app.simulations[segment].N &&
        segment_snapshot.nvars == app.sim.model.nvars &&
        length(app.plot_panel.segment_observables[segment]) ==
            segment_snapshot.nvars
        for (segment, segment_snapshot) in enumerate(snapshot.segments)
    )
end


function refresh_app_observables_from_snapshot!(
    app::AppState,
    snapshot::PartitionSnapshot,
)
    app.time_obs[] = snapshot.t
    app.dt_obs[] = snapshot.dt
    app.step_counter_obs[] = snapshot.steps

    return nothing
end


function refresh_app_from_snapshot!(
    app::AppState,
    snapshot::PartitionSnapshot,
)
    snapshot_matches_current_app(app, snapshot) ||
        return nothing

    refresh_plot_panel_from_snapshots!(app.plot_panel, snapshot.segments)
    refresh_app_observables_from_snapshot!(app, snapshot)

    return nothing
end


function refresh_app_observables_from_live_state!(app::AppState)
    app.time_obs[] = current_display_time(app.sim)
    app.dt_obs[] = current_internal_dt(app.sim)
    app.dtmax_obs[] = current_dtmax(app.sim)
    app.step_counter_obs[] = app.sim.step_counter[]

    return nothing
end


function refresh_app_from_live_state!(app::AppState)
    refresh_plot_panel!(app.plot_panel, app.simulations)
    refresh_app_observables_from_live_state!(app)

    return nothing
end


# ============================================================
# Worker control
# ============================================================

function request_stop_worker!(app::AppState)
    app.worker_running[] = false
    app.running[] = false
    app.synchronization_running[] = false

    for runtime in app.segment_runtimes
        runtime.running[] = false
    end

    return nothing
end


function wait_for_worker!(app::AppState)
    for runtime in app.segment_runtimes
        task = runtime.task_ref[]

        if task !== nothing && task !== current_task() && !istaskdone(task)
            wait(task)
        end
    end

    sync_task = app.synchronization_task_ref[]

    if sync_task !== nothing &&
       sync_task !== current_task() &&
       !istaskdone(sync_task)
        wait(sync_task)
    end

    return nothing
end


function stop_worker!(app::AppState; wait::Bool = false)
    request_stop_worker!(app)

    if wait
        wait_for_worker!(app)
    end

    return nothing
end


function refresh_ui_from_latest_snapshots!(app::AppState)
    snapshots = SimulationSnapshot[]

    for runtime in app.segment_runtimes
        lock(runtime.snapshot_lock)

        try
            runtime.latest_snapshot[] === nothing ||
                push!(snapshots, runtime.latest_snapshot[])
        finally
            unlock(runtime.snapshot_lock)
        end
    end

    if length(snapshots) == length(app.simulations) && !isempty(snapshots)
        try
            snapshot = partition_snapshot_from_segments(
                snapshots,
                app.generation[],
            )
            refresh_app_from_snapshot!(app, snapshot)
        catch err
            @error "Error while refreshing UI from snapshot." exception = (err, catch_backtrace())
        end
    end

    workers_active = any(runtime.running[] for runtime in app.segment_runtimes)
    app.worker_running[] = workers_active
    app.running[] = workers_active

    if !app.synchronization_running[] && length(snapshots) > 1
        times = [snapshot.t for snapshot in snapshots]
        tolerance = max(maximum(abs, times), 1.0) * 1e-10
        app.synchronization_status[] =
            maximum(times) - minimum(times) <= tolerance ?
            "Synchronized" :
            "Synchronize"
    elseif length(snapshots) <= 1 && !app.synchronization_running[]
        app.synchronization_status[] = "Synchronized"
    end

    return nothing
end


function start_ui_snapshot_poller!(
    app::AppState;
    refresh_interval::Float64 = 1 / 30,
)
    if app.ui_task_ref[] !== nothing && !istaskdone(app.ui_task_ref[])
        return nothing
    end

    app.ui_task_ref[] = @async begin
        while true
            refresh_ui_from_latest_snapshots!(app)
            yield()
            sleep(refresh_interval)
        end
    end

    return nothing
end


function start_worker!(
    app::AppState;
    steps_per_frame::Int = 5,
    sleep_time::Float64 = 0.001,
)
    app.synchronization_running[] && return nothing

    # When the simulation starts, hide perturbation previews.
    clear_perturbation_previews!(app.plot_panel)

    app.worker_running[] = true
    app.running[] = true

    for segment in eachindex(app.simulations)
        sim = app.simulations[segment]
        runtime = app.segment_runtimes[segment]
        existing_task = runtime.task_ref[]

        if existing_task !== nothing && !istaskdone(existing_task)
            runtime.running[] = true
            continue
        end

        runtime.running[] = true
        generation = app.generation[]
        runtime.task_ref[] = Threads.@spawn begin
            while runtime.running[]
                snapshot = nothing
                lock(runtime.lock)

                try
                    step_simulation!(sim, steps_per_frame)
                    snapshot = make_snapshot(sim, generation)
                catch err
                    @error "Critical error in segment worker." segment exception = (err, catch_backtrace())
                    runtime.running[] = false
                finally
                    unlock(runtime.lock)
                end

                if snapshot !== nothing
                    lock(runtime.snapshot_lock)

                    try
                        runtime.latest_snapshot[] = snapshot
                    finally
                        unlock(runtime.snapshot_lock)
                    end
                end

                yield()
                sleep(sleep_time)
            end
        end
    end

    active_tasks = [
        runtime.task_ref[]
        for runtime in app.segment_runtimes
        if runtime.task_ref[] !== nothing
    ]
    app.worker_task_ref[] = isempty(active_tasks) ? nothing : first(active_tasks)

    return nothing
end


function stop_segment_workers!(app::AppState; wait::Bool = false)
    app.worker_running[] = false
    app.running[] = false

    for runtime in app.segment_runtimes
        runtime.running[] = false
    end

    if wait
        for runtime in app.segment_runtimes
            task = runtime.task_ref[]

            if task !== nothing && task !== current_task() && !istaskdone(task)
                Base.wait(task)
            end
        end
    end

    return nothing
end


# ============================================================
# Controlled simulation mutations
# ============================================================

function make_locked_snapshot(app::AppState)
    # The caller should already hold app.simlock.

    return make_partition_snapshot(app.simulations, app.generation[])
end


function store_runtime_snapshots!(
    app::AppState,
    snapshots::Vector{SimulationSnapshot},
)
    length(snapshots) == length(app.segment_runtimes) ||
        error("Snapshot count does not match segment runtime count.")

    for segment in eachindex(snapshots)
        runtime = app.segment_runtimes[segment]
        lock(runtime.snapshot_lock)

        try
            runtime.latest_snapshot[] = snapshots[segment]
        finally
            unlock(runtime.snapshot_lock)
        end
    end

    return nothing
end


function step_once_app!(app::AppState)
    # Advance the simulation by one solver step.
    #
    # This is currently not attached to a button, but keeping it is useful
    # for debugging.

    if app.worker_running[]
        return nothing
    end

    snapshots = SimulationSnapshot[]

    for segment in eachindex(app.simulations)
        runtime = app.segment_runtimes[segment]
        lock(runtime.lock)

        try
            step_simulation!(app.simulations[segment])
            push!(snapshots, make_snapshot(app.simulations[segment], app.generation[]))
        finally
            unlock(runtime.lock)
        end
    end

    refresh_app_from_snapshot!(
        app,
        partition_snapshot_from_segments(snapshots, app.generation[]),
    )
    store_runtime_snapshots!(app, snapshots)

    return nothing
end


function set_dtmax_app!(app::AppState, new_dtmax::Float64)
    if app.synchronization_running[]
        app.dtmax_obs[] = new_dtmax
        return nothing
    end

    snapshots = SimulationSnapshot[]

    for segment in eachindex(app.simulations)
        runtime = app.segment_runtimes[segment]
        lock(runtime.lock)

        try
            set_dtmax!(app.simulations[segment], new_dtmax)
            snapshot = make_snapshot(app.simulations[segment], app.generation[])
            push!(snapshots, snapshot)

            lock(runtime.snapshot_lock)
            try
                runtime.latest_snapshot[] = snapshot
            finally
                unlock(runtime.snapshot_lock)
            end
        finally
            unlock(runtime.lock)
        end
    end

    app.dtmax_obs[] = new_dtmax
    refresh_app_from_snapshot!(
        app,
        partition_snapshot_from_segments(snapshots, app.generation[]),
    )

    return nothing
end


function with_worker_paused!(
    app::AppState,
    f::Function;
    restart_if_was_running::Bool = true,
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
)
    was_running = app.worker_running[]

    if was_running
        stop_worker!(app; wait = true)
    end

    snapshot = nothing

    lock(app.simlock)

    try
        f()
        snapshot = make_locked_snapshot(app)

    finally
        unlock(app.simlock)
    end

    clear_snapshot_buffer!(app.snapshot_buffer)

    if snapshot !== nothing
        store_runtime_snapshots!(app, snapshot.segments)
        refresh_app_from_snapshot!(app, snapshot)
    end

    if was_running && restart_if_was_running
        start_worker!(
            app;
            steps_per_frame = steps_per_frame,
            sleep_time = worker_sleep_time,
        )
    end

    return nothing
end


# ============================================================
# Reset and constant initial conditions
# ============================================================

function reset_initial_condition!(sim::SimulationState)
    # Reset the solution to the model-defined initial condition.

    U0 = zeros(Float64, sim.N, sim.model.nvars)

    sim.model.initialize!(
        U0,
        sim.x,
        sim.params,
    )

    ynew = copy(vec(U0))

    restart_after_manual_change!(sim, ynew)

    sim.time_offset[] = 0.0
    sim.step_counter[] = 0

    return nothing
end


function reset_initial_condition_app!(
    app::AppState;
    plot_grid::GridLayout,
    title_obs,
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
)
    stop_worker!(app; wait = true)
    domain_length_scale = app.plot_panel.domain_length_scale

    lock(app.simlock)

    try
        app.generation[] += 1
        clear_snapshot_buffer!(app.snapshot_buffer)

        app.sim = create_simulation_state(
            app.sim.model;
            N = app.initial_N,
            boundary_condition = app.initial_boundary_condition,
            dtmax = current_dtmax(app.sim),
        )
        set_diffusion_scale!(app.sim, domain_length_scale^2)
        app.simulations = SimulationState[app.sim]
        app.segment_runtimes = SegmentRuntime[empty_segment_runtime()]

        rebuild_plot_panel_for_partition!(
            app,
            plot_grid;
            title_obs = title_obs,
            domain_length_scale = domain_length_scale,
        )

        app.time_obs[] = 0.0
        app.step_counter_obs[] = 0
    finally
        unlock(app.simlock)
    end

    return nothing
end


# ============================================================
# In-memory save / restore
# ============================================================

function save_simulation_state!(app::AppState)
    length(app.simulations) == length(app.segment_runtimes) ||
        error("Segment runtime count does not match the simulation partition.")

    saved_segments = SavedSegmentState[]

    # Holding simlock freezes only topology-changing UI operations. Each live
    # solver keeps running until its own short runtime lock is acquired for a
    # consistent copy of that segment.
    lock(app.simlock)

    try
        generation = app.generation[]
        model_id = app.sim.model.id

        for segment in eachindex(app.simulations)
            sim = app.simulations[segment]
            runtime = app.segment_runtimes[segment]
            lock(runtime.lock)

            try
                app.generation[] == generation ||
                    error("Simulation topology changed while saving state.")
                sim.model.id == model_id ||
                    error("All saved segments must use the same model.")

                dtmax = current_dtmax(sim)
                if !isfinite(dtmax) || dtmax <= 0
                    dtmax = app.dtmax_obs[]
                end

                push!(
                    saved_segments,
                    SavedSegmentState(
                        model_id,
                        copy(sim.x),
                        sim.dx,
                        copy(sim.integrator_ref[].u),
                        deepcopy(sim.params),
                        sim.boundary_condition,
                        current_display_time(sim),
                        current_internal_dt(sim),
                        dtmax,
                        sim.step_counter[],
                    ),
                )
            finally
                unlock(runtime.lock)
            end
        end

        app.saved_state[] = SavedSimulationState(
            model_id,
            saved_segments,
            app.plot_panel.domain_length_scale,
            app.dtmax_obs[],
            app.initial_N,
            app.initial_boundary_condition,
        )
    finally
        unlock(app.simlock)
    end

    return app.saved_state[]
end


function clear_saved_simulation_state!(app::AppState)
    app.saved_state[] = nothing
    return nothing
end


function restore_saved_simulation_state_app!(
    app::AppState;
    plot_grid::GridLayout,
    title_obs,
    reltol::Float64 = 1e-5,
    abstol::Float64 = 1e-7,
)
    saved = app.saved_state[]
    saved === nothing && error("No saved simulation state is available.")
    isempty(saved.segments) && error("Saved simulation state has no segments.")
    saved.model_id == app.sim.model.id ||
        error("The saved state belongs to a different model.")

    stop_worker!(app; wait = true)
    restored_simulations = SimulationState[]

    lock(app.simlock)

    try
        app.generation[] += 1
        clear_snapshot_buffer!(app.snapshot_buffer)
        model = app.sim.model

        for segment in saved.segments
            segment.model_id == model.id ||
                error("Saved segment belongs to a different model.")
            isfinite(segment.dtmax) && segment.dtmax > 0 ||
                error("Saved segment has an invalid maximum time step.")

            sim = create_simulation_state_from_data(
                model,
                copy(segment.x),
                segment.dx,
                copy(segment.y),
                deepcopy(segment.params);
                boundary_condition = segment.boundary_condition,
                displayed_time = segment.t,
                dtmax = segment.dtmax,
                reltol = reltol,
                abstol = abstol,
            )
            sim.step_counter[] = segment.steps

            if isfinite(segment.dt) && segment.dt > 0
                set_proposed_dt!(
                    sim.integrator_ref[],
                    min(segment.dt, segment.dtmax),
                )
            end

            push!(restored_simulations, sim)
        end

        app.simulations = restored_simulations
        app.sim = first(restored_simulations)
        app.segment_runtimes = [
            empty_segment_runtime() for _ in restored_simulations
        ]
        app.initial_N = saved.initial_N
        app.initial_boundary_condition = saved.initial_boundary_condition
        app.dtmax_obs[] = saved.requested_dtmax
        app.running[] = false
        app.worker_running[] = false
        app.synchronization_running[] = false

        rebuild_plot_panel_for_partition!(
            app,
            plot_grid;
            title_obs = title_obs,
            domain_length_scale = saved.domain_length_scale,
        )

        snapshot = make_partition_snapshot(
            app.simulations,
            app.generation[],
        )
        store_runtime_snapshots!(app, snapshot.segments)
        refresh_app_from_snapshot!(app, snapshot)

        times = [segment.t for segment in saved.segments]
        tolerance = max(maximum(abs, times), 1.0) * 1e-10
        app.synchronization_status[] =
            maximum(times) - minimum(times) <= tolerance ?
            "Synchronized" :
            "Synchronize"
    finally
        unlock(app.simlock)
    end

    return nothing
end


function catch_up_segment_to_time!(
    app::AppState,
    segment::Int,
    target_time::Float64,
    ;
    allow_cancel::Bool,
)
    sim = app.simulations[segment]
    runtime = app.segment_runtimes[segment]
    tolerance = max(abs(target_time), 1.0) * 1e-11

    lock(runtime.lock)

    try
        remaining = target_time - current_display_time(sim)

        if remaining > tolerance
            integrator = sim.integrator_ref[]
            integrator.opts.dtmax = Inf
            add_tstop!(integrator, integrator.t + remaining)

            while current_display_time(sim) < target_time - tolerance
                allow_cancel && !app.synchronization_running[] && break
                step_simulation!(sim)

                if sim.step_counter[] % 5 == 0
                    snapshot = make_snapshot(sim, app.generation[])
                    lock(runtime.snapshot_lock)
                    try
                        runtime.latest_snapshot[] = snapshot
                    finally
                        unlock(runtime.snapshot_lock)
                    end
                end
            end
        end
    finally
        set_dtmax!(sim, app.dtmax_obs[])
        snapshot = make_snapshot(sim, app.generation[])
        lock(runtime.snapshot_lock)
        try
            runtime.latest_snapshot[] = snapshot
        finally
            unlock(runtime.snapshot_lock)
        end
        unlock(runtime.lock)
    end

    return nothing
end


function synchronize_segment_indices_blocking!(
    app::AppState,
    indices::AbstractVector{Int};
    allow_cancel::Bool,
)
    isempty(indices) && return nothing
    all(index -> 1 <= index <= length(app.simulations), indices) || return nothing
    target_time = maximum(current_display_time(app.simulations[index]) for index in indices)
    tasks = Task[]

    for index in indices
        if current_display_time(app.simulations[index]) < target_time
            task = Threads.@spawn begin
                catch_up_segment_to_time!(
                    app,
                    index,
                    target_time;
                    allow_cancel = allow_cancel,
                )
            end
            push!(tasks, task)
        else
            set_dtmax!(app.simulations[index], app.dtmax_obs[])
        end
    end

    wait.(tasks)

    for index in indices
        set_dtmax!(app.simulations[index], app.dtmax_obs[])
    end

    return nothing
end


function synchronize_domains_app!(app::AppState)
    length(app.simulations) <= 1 && begin
        app.synchronization_status[] = "Synchronized"
        return nothing
    end

    existing_task = app.synchronization_task_ref[]
    existing_task !== nothing && !istaskdone(existing_task) && return nothing

    app.synchronization_running[] = true
    app.synchronization_status[] = "Synchronizing..."

    app.synchronization_task_ref[] = @async begin
        stop_segment_workers!(app; wait = true)

        try
            synchronize_segment_indices_blocking!(
                app,
                collect(eachindex(app.simulations));
                allow_cancel = true,
            )

            if app.synchronization_running[]
                snapshot = make_partition_snapshot(app.simulations, app.generation[])
                store_runtime_snapshots!(app, snapshot.segments)
                refresh_app_from_snapshot!(app, snapshot)
                app.synchronization_status[] = "Synchronized"
            end
        catch err
            @error "Error while synchronizing domain solvers." exception = (err, catch_backtrace())
            app.synchronization_status[] = "Synchronize"
        finally
            app.synchronization_running[] = false
        end
    end

    return nothing
end


function set_constant_initial_condition!(
    sim::SimulationState,
    values::AbstractVector{<:Real},
)
    # Set every variable to a spatially constant value.

    length(values) == sim.model.nvars ||
        error("Expected $(sim.model.nvars) values, got $(length(values)).")

    ynew = copy(sim.integrator_ref[].u)
    U = reshape(ynew, sim.N, sim.model.nvars)

    for j in 1:sim.model.nvars
        U[:, j] .= Float64(values[j])
    end

    restart_after_manual_change!(sim, ynew)

    sim.step_counter[] = 0

    return nothing
end


function set_constant_initial_condition_app!(
    app::AppState,
    values::AbstractVector{<:Real};
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
)
    with_worker_paused!(
        app,
        () -> begin
            for sim in app.simulations
                set_constant_initial_condition!(sim, values)
            end
        end;
        restart_if_was_running = true,
        steps_per_frame = steps_per_frame,
        worker_sleep_time = worker_sleep_time,
    )

    return nothing
end


function set_single_constant_initial_condition!(
    sim::SimulationState;
    variable::Int,
    value::Real,
)
    variable >= 1 && variable <= sim.model.nvars ||
        error("Invalid variable index: $(variable).")

    ynew = copy(sim.integrator_ref[].u)
    U = reshape(ynew, sim.N, sim.model.nvars)

    U[:, variable] .= Float64(value)

    restart_after_manual_change!(sim, ynew)

    sim.step_counter[] = 0

    return nothing
end


function set_single_constant_initial_condition_app!(
    app::AppState;
    variable::Int,
    value::Real,
    segment::Union{Nothing, Int} = nothing,
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
)
    if segment !== nothing
        1 <= segment <= length(app.simulations) ||
            error("Invalid simulation segment: $segment.")
    end

    with_worker_paused!(
        app,
        () -> begin
            target_segments = segment === nothing ?
                eachindex(app.simulations) :
                segment:segment

            for segment_index in target_segments
                set_single_constant_initial_condition!(
                    app.simulations[segment_index];
                    variable = variable,
                    value = value,
                )
            end
        end;
        restart_if_was_running = true,
        steps_per_frame = steps_per_frame,
        worker_sleep_time = worker_sleep_time,
    )

    return nothing
end


# ============================================================
# Active spatial profile switching
# ============================================================

function set_active_spatial_profile_set!(
    sim::SimulationState,
    index::Int,
)
    # Change the active spatial profile set used by the model RHS.
    #
    # The active profile set index is stored as a hidden parameter.
    # The model reaction can then use active spatial profiles as p.ρ, p.source, etc.

    nsets = length(sim.model.spatial_profile_sets)

    nsets > 0 ||
        return nothing

    index >= 1 && index <= nsets ||
        error("Invalid spatial profile set index: $(index).")

    sim.params[ACTIVE_SPATIAL_PROFILE_SET_PARAM] = Float64(index)

    # The RHS has changed, so restart the integrator with the current state.
    ynew = copy(sim.integrator_ref[].u)

    restart_after_manual_change!(sim, ynew)

    return nothing
end


function set_active_spatial_profile_set_app!(
    app::AppState,
    index::Int;
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
)
    # UI-safe wrapper for changing the active spatial profile set.
    #
    # If the simulation is running, pause it, change the profile, restart the
    # integrator, refresh the UI, and resume the worker.

    with_worker_paused!(
        app,
        () -> begin
            for sim in app.simulations
                nsets = length(sim.model.spatial_profile_sets)
                1 <= index <= nsets ||
                    error("Invalid spatial profile set index: $(index).")
                sim.params[ACTIVE_SPATIAL_PROFILE_SET_PARAM] = Float64(index)
            end

            refresh_partition_spatial_profile_overrides!(app.simulations)

            for sim in app.simulations
                restart_after_manual_change!(
                    sim,
                    copy(sim.integrator_ref[].u),
                )
            end
        end;
        restart_if_was_running = true,
        steps_per_frame = steps_per_frame,
        worker_sleep_time = worker_sleep_time,
    )

    return nothing
end




# ============================================================
# Local perturbation application
# ============================================================

function apply_local_perturbation!(
    sim::SimulationState;
    variable::Int,
    center::Float64,
    width::Float64,
    height::Float64,
    random_mode::Bool,
)
    # Older direct local perturbation function.
    #
    # The current UI mainly uses apply_local_perturbation_increment!,
    # because it lets the preview and the applied random perturbation match.

    variable >= 1 && variable <= sim.model.nvars ||
        error("Invalid variable index: $(variable).")

    width > 0 ||
        error("Perturbation width must be positive.")

    ynew = copy(sim.integrator_ref[].u)
    U = reshape(ynew, sim.N, sim.model.nvars)

    mask = local_perturbation_mask(
        sim.x,
        center,
        width;
        boundary_condition = sim.boundary_condition,
    )

    if random_mode
        for i in eachindex(mask)
            if mask[i]
                U[i, variable] += height * rand()
            end
        end
    else
        U[mask, variable] .+= height
    end

    restart_after_manual_change!(sim, ynew)

    sim.step_counter[] = 0

    return nothing
end


function apply_local_perturbation_app!(
    app::AppState;
    variable::Int,
    center::Float64,
    width::Float64,
    height::Float64,
    random_mode::Bool,
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
)
    with_worker_paused!(
        app,
        () -> apply_local_perturbation!(
            app.sim;
            variable = variable,
            center = center,
            width = width,
            height = height,
            random_mode = random_mode,
        );
        restart_if_was_running = true,
        steps_per_frame = steps_per_frame,
        worker_sleep_time = worker_sleep_time,
    )

    clear_perturbation_preview!(app.plot_panel, 1, variable)

    return nothing
end


function apply_local_perturbation_increment!(
    sim::SimulationState;
    variable::Int,
    increment::AbstractVector{<:Real},
)
    # Apply a precomputed perturbation increment to one variable.
    #
    # This is used by the local perturbation preview system. In random mode,
    # the random vector is generated at preview time and then exactly the same
    # vector is applied here.

    variable >= 1 && variable <= sim.model.nvars ||
        error("Invalid variable index: $(variable).")

    length(increment) == sim.N ||
        error("Perturbation increment has wrong length.")

    ynew = copy(sim.integrator_ref[].u)
    U = reshape(ynew, sim.N, sim.model.nvars)

    U[:, variable] .+= Float64.(increment)

    restart_after_manual_change!(sim, ynew)

    sim.step_counter[] = 0

    return nothing
end


function apply_local_perturbation_increment_app!(
    app::AppState;
    segment::Int,
    variable::Int,
    increment::AbstractVector{<:Real},
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
)
    with_worker_paused!(
        app,
        () -> apply_local_perturbation_increment!(
            app.simulations[segment];
            variable = variable,
            increment = increment,
        );
        restart_if_was_running = true,
        steps_per_frame = steps_per_frame,
        worker_sleep_time = worker_sleep_time,
    )

    clear_perturbation_preview!(app.plot_panel, segment, variable)

    return nothing
end


# ============================================================
# Plot panel rebuilding
# ============================================================

function switch_model_app!(
    app::AppState,
    plot_grid::GridLayout,
    model::ModelSpec;
    N::Int,
    dtmax::Float64,
    reltol::Float64,
    abstol::Float64,
    boundary_condition::Symbol,
    title_obs,
    model_name_obs::Observable{String},
    bc_name_obs = nothing,
)
    # Switch to another already-loaded model.

    stop_worker!(app; wait = true)

    lock(app.simlock)

    try
        app.generation[] = app.generation[] + 1
        clear_snapshot_buffer!(app.snapshot_buffer)
        clear_saved_simulation_state!(app)

        selected_boundary_condition = isnothing(model.default_boundary_condition) ?
            boundary_condition : model.default_boundary_condition

        app.sim = create_simulation_state(
            model;
            N = N,
            dtmax = dtmax,
            reltol = reltol,
            abstol = abstol,
            boundary_condition = selected_boundary_condition,
        )
        app.simulations = SimulationState[app.sim]
        app.segment_runtimes = SegmentRuntime[empty_segment_runtime()]
        app.initial_N = N
        app.initial_boundary_condition = selected_boundary_condition

        clear_plot_panel!(app.plot_panel)
        reset_plot_grid_layout!(plot_grid)

        app.plot_panel = build_plot_panel!(
            plot_grid,
            app;
            title_obs = title_obs,
        )

        refresh_app_from_live_state!(app)

        # Important:
        #
        # Update model_name_obs only after app.sim and app.plot_panel
        # are already consistent.
        #
        # This observable triggers callbacks in ControlPanel.jl, including
        # diffusion rescaling. If it is triggered too early, the new model
        # may have 3 variables while the old plot panel still has 2 axes.
        model_name_obs[] = model.display_name

        if bc_name_obs !== nothing
            bc_name_obs[] = boundary_condition_label(selected_boundary_condition)
        end

    finally
        unlock(app.simlock)
    end

    return nothing
end

function switch_boundary_condition_app!(
    app::AppState,
    plot_grid::GridLayout,
    boundary_condition::Symbol;
    N::Int,
    dtmax::Float64,
    reltol::Float64,
    abstol::Float64,
    title_obs,
    bc_name_obs::Observable{String},
)
    # Switch boundary condition for the currently selected model.
    #
    # This rebuilds the grid, the Laplacian, the initial condition,
    # and the solver.

    validate_boundary_condition(boundary_condition)

    stop_worker!(app; wait = true)

    lock(app.simlock)

    try
        app.generation[] = app.generation[] + 1
        clear_snapshot_buffer!(app.snapshot_buffer)

        current_model = app.sim.model

        app.sim = create_simulation_state(
            current_model;
            N = N,
            dtmax = dtmax,
            reltol = reltol,
            abstol = abstol,
            boundary_condition = boundary_condition,
        )
        app.simulations = SimulationState[app.sim]
        app.segment_runtimes = SegmentRuntime[empty_segment_runtime()]
        app.initial_N = N
        app.initial_boundary_condition = boundary_condition

        bc_name_obs[] = boundary_condition_label(boundary_condition)

        clear_plot_panel!(app.plot_panel)
        reset_plot_grid_layout!(plot_grid)

        app.plot_panel = build_plot_panel!(
            plot_grid,
            app;
            title_obs = title_obs,
        )

        refresh_app_from_live_state!(app)

    finally
        unlock(app.simlock)
    end

    return nothing
end

    # ============================================================
    # Diffusion rescaling
    # ============================================================

    function diffusion_parameter_names(model::ModelSpec)
        names = Symbol[]

        # Preferred convention:
        #
        #     variable u -> Du
        #     variable v -> Dv
        #     variable w -> Dw
        #
        for varname in model.varnames
            key = Symbol("D", varname)

            if haskey(model.default_params, key)
                push!(names, key)
            end
        end

        # Fallback: all parameters whose name starts with D.
        if isempty(names)
            for key in keys(model.default_params)
                if startswith(String(key), "D")
                    push!(names, key)
                end
            end
        end

        return sort(unique(names); by = String)
    end


    function diffusion_scale_label_string(sim::SimulationState, scale::Float64)
        names = diffusion_parameter_names(sim.model)

        if isempty(names)
            return "diffusion scale = $(@sprintf("%.2g", scale)) | no diffusion parameters found"
        end

        parts = String[]

        for name in names
            value = sim.params[name]
            push!(parts, "$(name) = $(@sprintf("%.2e", value))")
        end

        return "Length rescaled: $(@sprintf("%.2g", sqrt(scale))) "
    end


    function set_diffusion_scale!(sim::SimulationState, scale::Real)
        scale_float = Float64(scale)

        isfinite(scale_float) || error("Diffusion scale must be finite.")
        scale_float > 0.0 || error("Diffusion scale must be positive.")

        names = diffusion_parameter_names(sim.model)

        for name in names
            base_value = sim.model.default_params[name]
            sim.params[name] = base_value / scale_float
        end

        if haskey(sim.params, :domain_scale)
            sim.params[:domain_scale] = sqrt(scale_float)
        end

        # Restart the integrator from the current solution, with changed parameters.
        ynew = copy(sim.integrator_ref[].u)

        restart_after_manual_change!(
            sim,
            ynew,
        )

        return nothing
    end


    function set_diffusion_scale_app!(
        app::AppState,
        scale::Real;
        steps_per_frame::Int,
        worker_sleep_time::Float64,
    )
        with_worker_paused!(
            app,
            () -> begin
                for sim in app.simulations
                    set_diffusion_scale!(sim, scale)
                end
                set_plot_domain_scale!(app.plot_panel, app, scale)
            end;
            restart_if_was_running = true,
            steps_per_frame = steps_per_frame,
            worker_sleep_time = worker_sleep_time,
        )

        return nothing
    end
