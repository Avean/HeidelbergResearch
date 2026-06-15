# src/UI.jl

# ============================================================
# Main UI composition
# ============================================================
#
# The window is composed of three visual areas:
#
#     1. Top menu
#     2. Plot panel
#     3. Control panel
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
    N::Int = 1000,
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
    # Main window layout
    # --------------------------------------------------------

    fig = Figure(size = (1300, 820))

    top_menu_grid = GridLayout()
    plot_grid = GridLayout()
    control_grid = GridLayout()

    fig[1, 1:2] = top_menu_grid
    fig[2, 1] = plot_grid
    fig[2, 2] = control_grid

    rowsize!(fig.layout, 1, Fixed(50))
    colsize!(fig.layout, 1, Relative(0.68))
    colsize!(fig.layout, 2, Relative(0.32))

    # --------------------------------------------------------
    # Plot panel
    # --------------------------------------------------------

    plot_panel = empty_plot_panel()

    # --------------------------------------------------------
    # Application state
    # --------------------------------------------------------

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

    app.plot_panel = build_plot_panel!(
        plot_grid,
        app;
        title_obs = title_obs,
    )

    # --------------------------------------------------------
    # Top menu
    # --------------------------------------------------------

    build_top_menu!(
        top_menu_grid,
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
        dtmax0 = dtmax0,
        steps_per_frame = steps_per_frame,
        worker_sleep_time = worker_sleep_time,
        model_name_obs = model_name_obs,
    )

    # --------------------------------------------------------
    # Start
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