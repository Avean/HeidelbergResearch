# src/UI.jl

# ============================================================
# User interface
# ============================================================
#
# This file builds the GLMakie application window.
#
# Threaded architecture:
#
#     Worker thread:
#         - advances the DifferentialEquations.jl integrator,
#         - creates SimulationSnapshot objects,
#         - writes only to SnapshotBuffer.
#
#     UI task:
#         - reads the newest SimulationSnapshot,
#         - updates Makie Observables,
#         - updates axes and labels.
#
# The worker thread must never update GLMakie objects directly.
#
# ============================================================


# ============================================================
# Snapshot-based UI refresh
# ============================================================

function snapshot_matches_current_app(
    app::AppState,
    snapshot::SimulationSnapshot,
)
    # Check whether a snapshot still belongs to the currently displayed model.
    #
    # This prevents drawing an old snapshot after switching models.

    return snapshot.generation == app.generation[] &&
           snapshot.model_id == app.sim.model.id &&
           snapshot.N == app.sim.N &&
           snapshot.nvars == app.sim.model.nvars
end


function refresh_app_observables_from_snapshot!(
    app::AppState,
    snapshot::SimulationSnapshot,
)
    # Refresh scalar observables displayed in the UI using a snapshot.

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
    # Refresh plots and scalar UI observables from a thread-safe snapshot.
    #
    # This function should be called only from the UI side.

    snapshot_matches_current_app(app, snapshot) ||
        return nothing

    refresh_plot_panel_from_snapshot!(app.plot_panel, snapshot)
    refresh_app_observables_from_snapshot!(app, snapshot)

    return nothing
end


function refresh_app_observables_from_live_state!(app::AppState)
    # Refresh scalar observables from the live simulation state.
    #
    # Use this only when the worker is stopped or when the caller knows that
    # access to the simulation state is safe.

    app.time_obs[] = current_display_time(app.sim)
    app.dt_obs[] = current_internal_dt(app.sim)
    app.dtmax_obs[] = current_dtmax(app.sim)
    app.step_counter_obs[] = app.sim.step_counter[]

    return nothing
end


function refresh_app_from_live_state!(app::AppState)
    # Refresh plots and scalar UI observables directly from the live state.
    #
    # This is used at initialization and immediately after controlled manual
    # operations, not during threaded stepping.

    refresh_plot_panel!(app.plot_panel, app.sim)
    refresh_app_observables_from_live_state!(app)

    return nothing
end


# ============================================================
# Worker control
# ============================================================

function request_stop_worker!(app::AppState)
    # Ask the worker thread to stop.
    #
    # This does not wait until the worker actually finishes.

    app.worker_running[] = false
    app.running[] = false

    return nothing
end


function wait_for_worker!(app::AppState)
    # Wait until the worker task finishes.
    #
    # Important:
    # Do not call this while holding app.simlock.

    task = app.worker_task_ref[]

    if task !== nothing && !istaskdone(task)
        wait(task)
    end

    return nothing
end


function stop_worker!(app::AppState; wait::Bool = false)
    # Stop the worker.
    #
    # For buttons, wait = false is more responsive.
    # For model switching or solver mutation, wait = true is safer.

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
    # Start the UI-side task that polls the latest snapshot and updates Makie.
    #
    # This task runs on the UI side. It does not step the solver.

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
    # Start the threaded simulation worker.
    #
    # The worker advances the solver and publishes snapshots.
    # It must not touch GLMakie objects.

    if app.worker_task_ref[] !== nothing && !istaskdone(app.worker_task_ref[])
        app.worker_running[] = true
        app.running[] = true
        return nothing
    end

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
    # Create a snapshot while holding the simulation lock.
    #
    # The caller must already hold app.simlock.

    return make_snapshot(app.sim, app.generation[])
end


function step_once_app!(app::AppState)
    # Advance the simulation by one solver step.
    #
    # If the worker is running, this button does nothing. Manual stepping is
    # intended for the stopped state.

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
    # Change dtmax through the UI.
    #
    # This can be done while the worker is running. The solver state is locked
    # briefly, and the UI is refreshed from a snapshot.

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
)
    # Pause the worker, apply a mutation to the live simulation state,
    # refresh from a snapshot, and optionally restart the worker.
    #
    # This is used for manual perturbations.

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
        )
    end

    return nothing
end


function kick_all_app!(
    app::AppState;
    steps_per_frame::Int = 5,
)
    # Add a constant perturbation to the first variable.

    with_worker_paused!(
        app,
        () -> kick_all!(app.sim; amount = 1.0, variable = 1);
        restart_if_was_running = true,
        steps_per_frame = steps_per_frame,
    )

    return nothing
end


function sinusoidal_kick_app!(
    app::AppState;
    steps_per_frame::Int = 5,
)
    # Add a sinusoidal perturbation to the first variable.

    with_worker_paused!(
        app,
        () -> sinusoidal_kick!(
            app.sim;
            amount = 1.0,
            mode = 1,
            variable = 1,
            shift = 0.0,
            clip_to_nonnegative = false,
        );
        restart_if_was_running = true,
        steps_per_frame = steps_per_frame,
    )

    return nothing
end


function sinusoidal_kick_clipped_app!(
    app::AppState;
    steps_per_frame::Int = 5,
)
    # Add a shifted sinusoidal perturbation and clip the first variable
    # to nonnegative values.

    with_worker_paused!(
        app,
        () -> sinusoidal_kick!(
            app.sim;
            amount = 1.0,
            mode = 1,
            variable = 1,
            shift = 1.5,
            clip_to_nonnegative = true,
        );
        restart_if_was_running = true,
        steps_per_frame = steps_per_frame,
    )

    return nothing
end


# ============================================================
# Plot panel rebuilding
# ============================================================

function clear_plot_panel!(panel::PlotPanel)
    # Delete old axes before rebuilding the plot panel.
    #
    # This is needed when the selected model changes, because different
    # models can have different numbers of variables.

    for ax in panel.axes
        try
            delete!(ax)
        catch
            try
                ax.visible = false
            catch
            end
        end
    end

    empty!(panel.axes)
    empty!(panel.observables)

    return nothing
end


function switch_model_app!(
    app::AppState,
    plot_grid::GridLayout,
    model::ModelSpec;
    N::Int,
    dtmax::Float64,
    reltol::Float64,
    abstol::Float64,
    title_obs,
    model_name_obs::Observable{String},
)
    # Switch to another already-loaded model.
    #
    # The worker is stopped, old snapshots are discarded, and a new generation
    # number is assigned so stale snapshots cannot be drawn.

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
        )

        model_name_obs[] = model.display_name

        clear_plot_panel!(app.plot_panel)

        app.plot_panel = build_plot_panel!(
            plot_grid,
            app.sim;
            title_obs = title_obs,
        )

        refresh_app_from_live_state!(app)

    finally
        unlock(app.simlock)
    end

    return nothing
end


# ============================================================
# Main UI
# ============================================================

function run_app(;
    N::Int = 1000,
    dtmax0::Float64 = 1e-2,
    reltol::Float64 = 1e-5,
    abstol::Float64 = 1e-7,
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
    ui_refresh_interval::Float64 = 1 / 30,
)
    # Start the full GLMakie application.

    GLMakie.activate!()

    if Threads.nthreads() == 1
        @warn """
        Julia is running with only one thread.

        The application will still work, but the simulation worker cannot run
        on a separate thread. Start Julia with for example:

            julia --threads=2 main.jl

        or set JULIA_NUM_THREADS before starting Julia.
        """
    end

    # --------------------------------------------------------
    # Models
    # --------------------------------------------------------

    registry = MODEL_REGISTRY
    labels = model_labels(registry)

    isempty(labels) &&
        error("No models found in MODEL_REGISTRY.")

    first_label = first(labels)
    first_model = get_model(registry, first_label)

    # --------------------------------------------------------
    # Initial simulation
    # --------------------------------------------------------

    sim = create_simulation_state(
        first_model;
        N = N,
        dtmax = dtmax0,
        reltol = reltol,
        abstol = abstol,
    )

    # --------------------------------------------------------
    # Observables
    # --------------------------------------------------------

    running_obs = Observable(false)
    dtmax_obs = Observable(dtmax0)
    dt_obs = Observable(current_internal_dt(sim))
    time_obs = Observable(current_display_time(sim))
    step_counter_obs = Observable(sim.step_counter[])
    model_name_obs = Observable(first_model.display_name)

    title_obs = lift(
        model_name_obs,
        time_obs,
        dtmax_obs,
        dt_obs,
        running_obs,
        step_counter_obs,
    ) do model_name, t, dtmax, dt, running, steps

        return "$(model_name) | TRBDF2 | running = $(running) | t = $(round(t; digits = 2)) | max dt = $(dtmax) | current dt = $(dt) | steps = $(steps)"
    end

    # --------------------------------------------------------
    # Figure layout
    # --------------------------------------------------------

    fig = Figure(size = (1200, 760))

    plot_grid = GridLayout()
    fig[1:18, 1] = plot_grid

    control_grid = GridLayout()
    fig[1:18, 2] = control_grid

    # --------------------------------------------------------
    # Initial plot panel
    # --------------------------------------------------------

    plot_panel = build_plot_panel!(
        plot_grid,
        sim;
        title_obs = title_obs,
    )

    app = AppState(
        sim,
        plot_panel,
        running_obs,
        dtmax_obs,
        dt_obs,
        time_obs,
        step_counter_obs,
        Threads.Atomic{Bool}(false),
        Ref{Union{Nothing, Task}}(nothing),
        Ref{Union{Nothing, Task}}(nothing),
        empty_snapshot_buffer(),
        Threads.Atomic{Int}(0),
        ReentrantLock(),
    )

    # --------------------------------------------------------
    # Model selector
    # --------------------------------------------------------

    Label(control_grid[1, 1], "Model", tellwidth = false)

    model_menu = Menu(
        control_grid[2, 1],
        options = labels,
        default = first_label,
        tellwidth = false,
    )

    on(model_menu.selection) do selected_label
        model = get_model(registry, selected_label)

        switch_model_app!(
            app,
            plot_grid,
            model;
            N = N,
            dtmax = current_dtmax(app.sim),
            reltol = reltol,
            abstol = abstol,
            title_obs = title_obs,
            model_name_obs = model_name_obs,
        )
    end

    # --------------------------------------------------------
    # Time-step controls
    # --------------------------------------------------------

    Label(control_grid[4, 1], "Choose max dt", tellwidth = false)

    dt_choices = [1e-5, 1e-2, 1.0, 1e5]

    for (i, dtchoice) in enumerate(dt_choices)
        b = Button(
            control_grid[4 + i, 1],
            label = "max dt = $(dtchoice)",
            tellwidth = false,
        )

        on(b.clicks) do _
            set_dtmax_app!(app, dtchoice)
        end
    end

    Label(control_grid[10, 1], "Custom max dt:", tellwidth = false)

    dtbox = Textbox(
        control_grid[11, 1],
        placeholder = "for example 2e-3",
        stored_string = string(dtmax0),
        tellwidth = false,
    )

    on(dtbox.stored_string) do s
        val = tryparse(Float64, s)

        if val !== nothing && isfinite(val) && val > 0
            set_dtmax_app!(app, val)
        end
    end

    # --------------------------------------------------------
    # Simulation buttons
    # --------------------------------------------------------

    bstart = Button(
        control_grid[13, 1],
        label = "Start",
        tellwidth = false,
    )

    bstop = Button(
        control_grid[14, 1],
        label = "Stop",
        tellwidth = false,
    )

    bone = Button(
        control_grid[15, 1],
        label = "One step",
        tellwidth = false,
    )

    on(bstart.clicks) do _
        start_worker!(
            app;
            steps_per_frame = steps_per_frame,
            sleep_time = worker_sleep_time,
        )
    end

    on(bstop.clicks) do _
        stop_worker!(app; wait = false)
    end

    on(bone.clicks) do _
        step_once_app!(app)
    end

    # --------------------------------------------------------
    # Perturbation buttons
    # --------------------------------------------------------

    bkick = Button(
        control_grid[17, 1],
        label = "kick: first variable += 1",
        tellwidth = false,
    )

    bsinkick = Button(
        control_grid[18, 1],
        label = "sin kick",
        tellwidth = false,
    )

    bsinkickclip = Button(
        control_grid[19, 1],
        label = "sin kick clipped",
        tellwidth = false,
    )

    on(bkick.clicks) do _
        kick_all_app!(
            app;
            steps_per_frame = steps_per_frame,
        )
    end

    on(bsinkick.clicks) do _
        sinusoidal_kick_app!(
            app;
            steps_per_frame = steps_per_frame,
        )
    end

    on(bsinkickclip.clicks) do _
        sinusoidal_kick_clipped_app!(
            app;
            steps_per_frame = steps_per_frame,
        )
    end

    # --------------------------------------------------------
    # Start UI poller and simulation worker
    # --------------------------------------------------------

    refresh_app_from_live_state!(app)

    display(fig)

    start_ui_snapshot_poller!(
        app;
        refresh_interval = ui_refresh_interval,
    )

    start_worker!(
        app;
        steps_per_frame = steps_per_frame,
        sleep_time = worker_sleep_time,
    )

    return app
end