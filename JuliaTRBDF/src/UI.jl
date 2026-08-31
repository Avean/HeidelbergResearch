# src/UI.jl

# ============================================================
# Main UI composition
# ============================================================
#
# The application is composed of two windows:
#
#     1. Main window: plot panel
#     2. Control window: model, boundary-condition, and simulation controls
#
# The actual logic is split into:
#
#     TopMenu.jl
#     PlotPanel.jl
#     ControlPanel.jl
#     UIRuntime.jl
#
# ============================================================


function run_app(;
    N::Int = 500,
    boundary_condition0::Symbol = :neumann,
    dtmax0::Float64 = 1e-2,
    reltol::Float64 = 1e-5,
    abstol::Float64 = 1e-7,
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
    ui_refresh_interval::Float64 = 1 / 30,
)
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

    validate_boundary_condition(boundary_condition0)

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
        boundary_condition = boundary_condition0,
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
    bc_name_obs = Observable(boundary_condition_label(boundary_condition0))

    title_obs = lift(
        model_name_obs,
        bc_name_obs,
        time_obs,
        dtmax_obs,
        dt_obs,
        running_obs,
        step_counter_obs,
    ) do model_name, bc_name, t, dtmax, dt, running, steps

        return "t = $(@sprintf("%.1e", t)) | steps = $(steps)"
    end

    # --------------------------------------------------------
    # Window layouts
    # --------------------------------------------------------

    main_fig = Figure(size = (1300, 820))
    control_fig = Figure(size = (620, 820))

    model_control_grid = GridLayout(tellheight = false)
    plot_grid = GridLayout(
        tellwidth = false,
        tellheight = false,
    )
    control_grid = GridLayout()

    main_fig[1, 1] = plot_grid
    control_fig[1, 1] = model_control_grid
    control_fig[2, 1] = control_grid

    rowsize!(main_fig.layout, 1, Auto(false, 1.0))
    rowsize!(control_fig.layout, 1, Fixed(50))
    rowsize!(control_fig.layout, 2, Auto(false, 1.0))
    rowgap!(main_fig.layout, 0)
    rowgap!(control_fig.layout, 0)
    colsize!(main_fig.layout, 1, Relative(1.0))
    colsize!(control_fig.layout, 1, Relative(1.0))

    # --------------------------------------------------------
    # Plot panel
    # --------------------------------------------------------

    plot_panel = empty_plot_panel()

    # --------------------------------------------------------
    # Application state
    # --------------------------------------------------------

    app = AppState(
        sim,
        SimulationState[sim],
        N,
        boundary_condition0,
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
        SegmentRuntime[empty_segment_runtime()],
        Threads.Atomic{Bool}(false),
        Ref{Union{Nothing, Task}}(nothing),
        Observable("Synchronized"),
        Ref{Union{Nothing, SavedSimulationState}}(nothing),
        true,
    )

    app.plot_panel = build_plot_panel!(
        plot_grid,
        app;
        title_obs = title_obs,
    )

    # --------------------------------------------------------
    # Model and boundary-condition controls
    # --------------------------------------------------------

    build_top_menu!(
        model_control_grid,
        app,
        plot_grid;
        registry = registry,
        labels = labels,
        first_label = first_label,
        N = N,
        reltol = reltol,
        abstol = abstol,
        title_obs = title_obs,
        model_name_obs = model_name_obs,
        bc_name_obs = bc_name_obs,
    )

    # --------------------------------------------------------
    # Control panel
    # --------------------------------------------------------

    build_control_panel!(
        control_grid,
        app;
        plot_grid = plot_grid,
        title_obs = title_obs,
        dtmax0 = dtmax0,
        steps_per_frame = steps_per_frame,
        worker_sleep_time = worker_sleep_time,
        model_name_obs = model_name_obs,
    )

    # --------------------------------------------------------
    # Start
    # --------------------------------------------------------

    refresh_app_from_live_state!(app)

    main_screen = display(
        GLMakie.Screen(title = "Reaction-Diffusion: plots"),
        main_fig,
    )

    control_screen = display(
        GLMakie.Screen(title = "Reaction-Diffusion: controls"),
        control_fig,
    )

    start_ui_snapshot_poller!(
        app;
        refresh_interval = ui_refresh_interval,
    )

    # --------------------------------------------------------
    # Windows closed
    # --------------------------------------------------------


    on(events(main_fig.scene).window_open) do is_open
        if !is_open
            if isopen(control_screen)
                close(control_screen)
            end

            task_before = app.worker_task_ref[]

            worker_was_active =
                app.worker_running[] ||
                (task_before !== nothing && !istaskdone(task_before))

            if worker_was_active
                @info "Window closed. Stopping simulation worker."
            else
                @info "Window closed. No active simulation worker."
            end

            stop_worker!(app; wait = true)

            task_after = app.worker_task_ref[]

            worker_stopped =
                !app.worker_running[] &&
                (task_after === nothing || istaskdone(task_after))

            if worker_stopped
                @info "Simulation worker stopped successfully."
            else
                @warn "Simulation worker may still be running."
            end
        end
    end

    return app
end



