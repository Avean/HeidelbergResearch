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
    snapshot::SimulationSnapshot,
)
    return snapshot.generation == app.generation[] &&
           snapshot.model_id == app.sim.model.id &&
           snapshot.N == app.sim.N &&
           snapshot.nvars == app.sim.model.nvars
end


function refresh_app_observables_from_snapshot!(
    app::AppState,
    snapshot::SimulationSnapshot,
)
    app.time_obs[] = snapshot.t
    app.dt_obs[] = snapshot.dt
    app.dtmax_obs[] = snapshot.dtmax
    app.step_counter_obs[] = snapshot.steps

    return nothing
end


function refresh_app_from_snapshot!(
    app::AppState,
    snapshot::SimulationSnapshot,
)
    snapshot_matches_current_app(app, snapshot) ||
        return nothing

    refresh_plot_panel_from_snapshot!(app.plot_panel, snapshot)
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
    refresh_plot_panel!(app.plot_panel, app.sim)
    refresh_app_observables_from_live_state!(app)

    return nothing
end


# ============================================================
# Worker control
# ============================================================

function request_stop_worker!(app::AppState)
    app.worker_running[] = false
    app.running[] = false

    return nothing
end


function wait_for_worker!(app::AppState)
    task = app.worker_task_ref[]

    if task !== nothing && !istaskdone(task)
        wait(task)
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


function start_ui_snapshot_poller!(
    app::AppState;
    refresh_interval::Float64 = 1 / 30,
)
    if app.ui_task_ref[] !== nothing && !istaskdone(app.ui_task_ref[])
        return nothing
    end

    app.ui_task_ref[] = @async begin
        while true
            snapshot = take_latest_snapshot!(app.snapshot_buffer)

            if snapshot !== nothing
                try
                    refresh_app_from_snapshot!(app, snapshot)
                catch err
                    @error "Error while refreshing UI from snapshot." exception = (err, catch_backtrace())
                end
            end

            # If the worker stopped because of an error, reflect that in the UI.
            if app.running[] &&
               !app.worker_running[] &&
               app.worker_task_ref[] !== nothing &&
               istaskdone(app.worker_task_ref[])
                app.running[] = false
            end

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
    if app.worker_task_ref[] !== nothing && !istaskdone(app.worker_task_ref[])
        app.worker_running[] = true
        app.running[] = true
        return nothing
    end

    # When the simulation starts, hide perturbation previews.
    clear_perturbation_previews!(app.plot_panel)

    app.worker_running[] = true
    app.running[] = true

    app.worker_task_ref[] = Threads.@spawn begin
        while app.worker_running[]
            snapshot = nothing

            lock(app.simlock)

            try
                step_simulation!(app.sim, steps_per_frame)

                snapshot = make_snapshot(
                    app.sim,
                    app.generation[],
                )

            catch err
                @error "Critical error in the simulation worker." exception = (err, catch_backtrace())
                app.worker_running[] = false

            finally
                unlock(app.simlock)
            end

            if snapshot !== nothing
                put_latest_snapshot!(app.snapshot_buffer, snapshot)
            end

            yield()
            sleep(sleep_time)
        end
    end

    return nothing
end


# ============================================================
# Controlled simulation mutations
# ============================================================

function make_locked_snapshot(app::AppState)
    # The caller should already hold app.simlock.

    return make_snapshot(app.sim, app.generation[])
end


function step_once_app!(app::AppState)
    # Advance the simulation by one solver step.
    #
    # This is currently not attached to a button, but keeping it is useful
    # for debugging.

    if app.worker_running[]
        return nothing
    end

    snapshot = nothing

    lock(app.simlock)

    try
        step_simulation!(app.sim)
        snapshot = make_locked_snapshot(app)

    finally
        unlock(app.simlock)
    end

    if snapshot !== nothing
        refresh_app_from_snapshot!(app, snapshot)
    end

    return nothing
end


function set_dtmax_app!(app::AppState, new_dtmax::Float64)
    snapshot = nothing

    lock(app.simlock)

    try
        set_dtmax!(app.sim, new_dtmax)
        snapshot = make_locked_snapshot(app)

    finally
        unlock(app.simlock)
    end

    if snapshot !== nothing
        refresh_app_from_snapshot!(app, snapshot)
    end

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

    sim.step_counter[] = 0

    return nothing
end


function reset_initial_condition_app!(
    app::AppState;
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
)
    with_worker_paused!(
        app,
        () -> reset_initial_condition!(app.sim);
        restart_if_was_running = true,
        steps_per_frame = steps_per_frame,
        worker_sleep_time = worker_sleep_time,
    )

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
        () -> set_constant_initial_condition!(app.sim, values);
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
        () -> set_active_spatial_profile_set!(
            app.sim,
            index,
        );
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

    clear_perturbation_preview!(app.plot_panel, variable)

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
    variable::Int,
    increment::AbstractVector{<:Real},
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
)
    with_worker_paused!(
        app,
        () -> apply_local_perturbation_increment!(
            app.sim;
            variable = variable,
            increment = increment,
        );
        restart_if_was_running = true,
        steps_per_frame = steps_per_frame,
        worker_sleep_time = worker_sleep_time,
    )

    clear_perturbation_preview!(app.plot_panel, variable)

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
)
    # Switch to another already-loaded model.

    stop_worker!(app; wait = true)

    lock(app.simlock)

    try
        app.generation[] = app.generation[] + 1
        clear_snapshot_buffer!(app.snapshot_buffer)

        app.sim = create_simulation_state(
            model;
            N = N,
            dtmax = dtmax,
            reltol = reltol,
            abstol = abstol,
            boundary_condition = boundary_condition,
        )

        model_name_obs[] = model.display_name

        clear_plot_panel!(app.plot_panel)

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

        bc_name_obs[] = boundary_condition_label(boundary_condition)

        clear_plot_panel!(app.plot_panel)

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