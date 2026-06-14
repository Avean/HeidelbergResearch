# src/UI.jl

# ============================================================
# User interface
# ============================================================
#
# This file builds the GLMakie application window.
#
# Responsibilities:
#
#     - create the Figure layout,
#     - create the model selection menu,
#     - create buttons and text boxes,
#     - connect UI callbacks to the simulation engine,
#     - run the asynchronous simulation loop,
#     - refresh the plot panel.
#
# Models are not loaded dynamically here.
# They are already loaded into MODEL_REGISTRY when the module starts.
#
# ============================================================


function refresh_app_observables!(app::AppState)
    # Refresh scalar observables displayed in the UI.

    app.time_obs[] = current_display_time(app.sim)
    app.dt_obs[] = current_internal_dt(app.sim)
    app.dtmax_obs[] = current_dtmax(app.sim)
    app.step_counter_obs[] = app.sim.step_counter[]

    return nothing
end


function refresh_app!(app::AppState)
    # Refresh both plots and scalar UI observables.

    refresh_plot_panel!(app.plot_panel, app.sim)
    refresh_app_observables!(app)

    return nothing
end


function stop_app!(app::AppState)
    # Stop the asynchronous simulation loop.

    app.running[] = false

    return nothing
end


function start_app!(
    app::AppState;
    steps_per_frame::Int = 5,
    sleep_time::Float64 = 0.001,
)
    # Start the asynchronous simulation loop.
    #
    # If the loop is already running, do nothing.

    if app.task_ref[] !== nothing && !istaskdone(app.task_ref[])
        app.running[] = true
        return nothing
    end

    app.running[] = true

    app.task_ref[] = @async begin
        while app.running[]
            lock(app.simlock)

            try
                step_simulation!(app.sim, steps_per_frame)
                refresh_app!(app)

            catch err
                @error "Critical error in the simulation loop." exception = (err, catch_backtrace())
                app.running[] = false

            finally
                unlock(app.simlock)
            end

            yield()
            sleep(sleep_time)
        end
    end

    return nothing
end


function step_once_app!(app::AppState)
    # Advance the simulation by one solver step.

    lock(app.simlock)

    try
        step_simulation!(app.sim)
        refresh_app!(app)

    finally
        unlock(app.simlock)
    end

    return nothing
end


function set_dtmax_app!(app::AppState, new_dtmax::Float64)
    # Change dtmax through the UI.

    lock(app.simlock)

    try
        set_dtmax!(app.sim, new_dtmax)
        refresh_app!(app)

    finally
        unlock(app.simlock)
    end

    return nothing
end


function kick_all_app!(app::AppState)
    # Add a constant perturbation to the first variable.

    lock(app.simlock)

    try
        kick_all!(app.sim; amount = 1.0, variable = 1)
        refresh_app!(app)

    finally
        unlock(app.simlock)
    end

    return nothing
end


function sinusoidal_kick_app!(app::AppState)
    # Add a sinusoidal perturbation to the first variable.

    lock(app.simlock)

    try
        sinusoidal_kick!(
            app.sim;
            amount = 1.0,
            mode = 1,
            variable = 1,
            shift = 0.0,
            clip_to_nonnegative = false,
        )

        refresh_app!(app)

    finally
        unlock(app.simlock)
    end

    return nothing
end


function sinusoidal_kick_clipped_app!(app::AppState)
    # Add a shifted sinusoidal perturbation and clip the first variable
    # to nonnegative values.

    lock(app.simlock)

    try
        sinusoidal_kick!(
            app.sim;
            amount = 1.0,
            mode = 1,
            variable = 1,
            shift = 1.5,
            clip_to_nonnegative = true,
        )

        refresh_app!(app)

    finally
        unlock(app.simlock)
    end

    return nothing
end


function clear_plot_panel!(panel::PlotPanel)
    # Delete old axes before rebuilding the plot panel.
    #
    # This is needed when the selected model changes, because different
    # models can have different numbers of variables.

    for ax in panel.axes
        try
            delete!(ax)
        catch
            # Fallback for Makie versions where delete!(ax) is unavailable.
            # Hiding is less clean, but prevents the application from crashing.
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
    # The model is taken from MODEL_REGISTRY, so no dynamic include happens here.

    lock(app.simlock)

    try
        stop_app!(app)

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

        refresh_app!(app)

    finally
        unlock(app.simlock)
    end

    return nothing
end


function run_app(;
    N::Int = 1000,
    dtmax0::Float64 = 1e-2,
    reltol::Float64 = 1e-5,
    abstol::Float64 = 1e-7,
    steps_per_frame::Int = 5,
)
    # Start the full GLMakie application.

    GLMakie.activate!()

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
        Ref{Union{Nothing, Task}}(nothing),
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
        start_app!(
            app;
            steps_per_frame = steps_per_frame,
        )
    end

    on(bstop.clicks) do _
        stop_app!(app)
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
        kick_all_app!(app)
    end

    on(bsinkick.clicks) do _
        sinusoidal_kick_app!(app)
    end

    on(bsinkickclip.clicks) do _
        sinusoidal_kick_clipped_app!(app)
    end

    # --------------------------------------------------------
    # Start
    # --------------------------------------------------------

    refresh_app!(app)

    display(fig)

    start_app!(
        app;
        steps_per_frame = steps_per_frame,
    )

    return app
end