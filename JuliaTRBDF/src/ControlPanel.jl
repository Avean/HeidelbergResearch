# src/ControlPanel.jl

# ============================================================
# Control panel
# ============================================================
#
# Right application area.
#
# Contains:
#
#     - dtmax slider,
#     - start/stop/one-step/reset buttons,
#     - constant initial condition section,
#     - perturbation buttons.
#
# ============================================================


function delete_control_items!(items)
    for item in items
        try
            delete!(item)
        catch
            try
                item.visible = false
            catch
            end
        end
    end

    empty!(items)

    return nothing
end


function rebuild_constant_initial_condition_panel!(
    grid::GridLayout,
    app::AppState,
    item_ref,
    textbox_ref;
    steps_per_frame::Int,
    worker_sleep_time::Float64,
)
    # Rebuild the constant initial condition section.
    #
    # This is needed because different models can have different variables.

    delete_control_items!(item_ref[])
    empty!(textbox_ref[])

    nvars = app.sim.model.nvars
    varnames = app.sim.model.varnames

    title = Label(
        grid[1, 1:nvars],
        "Constant initial condition",
        tellwidth = false,
    )

    push!(item_ref[], title)

    for j in 1:nvars
        variable_label = Label(
            grid[2, j],
            varnames[j],
            tellwidth = false,
        )

        textbox = Textbox(
            grid[3, j],
            placeholder = "0.0",
            stored_string = "0.0",
            tellwidth = false,
        )

        push!(item_ref[], variable_label)
        push!(item_ref[], textbox)
        push!(textbox_ref[], textbox)
    end

    apply_button = Button(
        grid[4, 1:nvars],
        label = "Apply constants",
        tellwidth = false,
    )

    push!(item_ref[], apply_button)

    on(apply_button.clicks) do _
        values = Float64[]

        for textbox in textbox_ref[]
            value = tryparse(Float64, textbox.stored_string[])

            if value === nothing
                @warn "Invalid constant initial condition value." textbox.stored_string[]
                return nothing
            end

            push!(values, value)
        end

        set_constant_initial_condition_app!(
            app,
            values;
            steps_per_frame = steps_per_frame,
            worker_sleep_time = worker_sleep_time,
        )
    end

    return nothing
end


function rebuild_latex_equation_panel!(
    grid::GridLayout,
    app::AppState,
    item_ref,
)
    delete_control_items!(item_ref[])

    equations = app.sim.model.latex_equations

    title = Label(
        grid[1, 2],
        "Equations",
        tellwidth = false,
        halign = :left,
    )

    push!(item_ref[], title)

    colsize!(grid, 1, Fixed(115))
    colsize!(grid, 2, Relative(1.0))
    colgap!(grid, 1)

    if isempty(equations)
        empty_label = Label(
            grid[2, 2],
            "No LaTeX equations specified for this model.",
            tellwidth = false,
            halign = :left,
        )

        push!(item_ref[], empty_label)

        return nothing
    end

    for (k, equation) in enumerate(equations)
        equation_label = Label(
            grid[k + 1, 2],
            latexstring(equation),
            tellwidth = false,
            halign = :left,
            fontsize = 22,
        )

        push!(item_ref[], equation_label)
    end

    return nothing
end

function build_control_panel!(
    grid::GridLayout,
    app::AppState;
    dtmax0::Float64,
    steps_per_frame::Int,
    worker_sleep_time::Float64,
    model_name_obs::Observable{String},
)
    # --------------------------------------------------------
    # Time-step controls
    # --------------------------------------------------------

    top_grid = GridLayout(
    tellheight = false,
    )

    grid[1, 1:2] = top_grid

    dt_info_label = Label(
        top_grid[1, 1],
        lift(app.dtmax_obs, app.dt_obs) do dtmax, dt
            return "max dt = $(@sprintf("%.1e", dtmax)) | current dt = $(@sprintf("%.1e", dt))"
        end,
        tellwidth = false,
    )

    dt_exponent0 = -2

    dt_slider = Slider(
        top_grid[2, 1],
        range = -5:1:5,
        startvalue = dt_exponent0,
        tellwidth = false,
    )

    on(dt_slider.value) do x
        x_int = Int(x)
        set_dtmax_app!(app, 10.0^x_int)
    end

    # --------------------------------------------------------
    # Main control row
    # --------------------------------------------------------




    colsize!(grid, 1, Relative(0.45))
    colsize!(grid, 2, Relative(0.55))

    # --------------------------------------------------------
    # Simulation buttons
    # --------------------------------------------------------

    label_simulation_control = Label(
        top_grid[3, 1],
        "Simulation controls",
        tellwidth = false,
    )
    
    action_grid = GridLayout(
        tellheight = false,
    )

    top_grid[4, 1] = action_grid

    bstart = Button(
        action_grid[1, 1],
        label = "Start",
        tellwidth = false,
    )

    bstop = Button(
        action_grid[1, 2],
        label = "Stop",
        tellwidth = false,
    )

    breset = Button(
        action_grid[1, 3],
        label = "Reset",
        tellwidth = false,
    )

    on(bstart.clicks) do _
        clear_perturbation_previews!(app.plot_panel)
        start_worker!(
            app;
            steps_per_frame = steps_per_frame,
            sleep_time = worker_sleep_time,
        )
    end

    on(bstop.clicks) do _
        stop_worker!(app; wait = true)

        update_all_perturbation_previews!(
            app;
            stop_simulation = false,
        )
    end


    on(breset.clicks) do _
        reset_initial_condition_app!(
            app;
            steps_per_frame = steps_per_frame,
            worker_sleep_time = worker_sleep_time,
        )
    end

    # --------------------------------------------------------
    # Latex Equation panel
    # --------------------------------------------------------

    equation_grid = GridLayout(
        tellheight = false,
    )

    grid[2, 1:2] = equation_grid

    equation_items = Ref(Any[])

    rebuild_latex_equation_panel!(
        equation_grid,
        app,
        equation_items,
    )


    # --------------------------------------------------------
    # Diffusion rescaling
    # --------------------------------------------------------

    diffusion_scale_obs = Observable(1.0)

    diffusion_scale_grid = GridLayout(
        tellheight = false,
    )

    grid[3, 1:2] = diffusion_scale_grid

    diffusion_scale_label_obs = Observable(
        diffusion_scale_label_string(app.sim, diffusion_scale_obs[])
    )

    diffusion_scale_label = Label(
        diffusion_scale_grid[1, 1],
        diffusion_scale_label_obs,
        tellwidth = false,
        halign = :center,
    )

    diffusion_scale_slider = Slider(
        diffusion_scale_grid[2, 1],
        range = 0.0:0.2:3.0,
        startvalue = 0.0,
        tellwidth = false,
    )

    function update_diffusion_scale_label!()
        diffusion_scale_label_obs[] =
            diffusion_scale_label_string(app.sim, diffusion_scale_obs[])

        return nothing
    end

    on(diffusion_scale_slider.value) do exponent
        scale = 10.0^Float64(exponent)

        diffusion_scale_obs[] = scale

        set_diffusion_scale_app!(
            app,
            scale;
            steps_per_frame = steps_per_frame,
            worker_sleep_time = worker_sleep_time,
        )

        update_diffusion_scale_label!()
    end

    # --------------------------------------------------------
    # Constant initial condition section
    # --------------------------------------------------------

    constant_ic_grid = GridLayout(
        tellheight = false,
    )

    grid[4, 1:2] = constant_ic_grid

    constant_ic_items = Ref(Any[])
    constant_ic_textboxes = Ref(Any[])

    rebuild_constant_initial_condition_panel!(
        constant_ic_grid,
        app,
        constant_ic_items,
        constant_ic_textboxes;
        steps_per_frame = steps_per_frame,
        worker_sleep_time = worker_sleep_time,
    )

    on(model_name_obs) do _
        rebuild_constant_initial_condition_panel!(
            constant_ic_grid,
            app,
            constant_ic_items,
            constant_ic_textboxes;
            steps_per_frame = steps_per_frame,
            worker_sleep_time = worker_sleep_time,
        )
        
        rebuild_latex_equation_panel!(
        equation_grid,
        app,
        equation_items,
        )

        set_diffusion_scale_app!(
        app,
        diffusion_scale_obs[];
        steps_per_frame = steps_per_frame,
        worker_sleep_time = worker_sleep_time,
        )

        update_diffusion_scale_label!()
    end


    return (;
        dt_info_label,
        label_simulation_control,
        dt_slider,
        action_grid,
        equation_grid,
        diffusion_scale_grid,
        diffusion_scale_slider,
        diffusion_scale_label,
        diffusion_scale_obs,
        constant_ic_grid,
        bstart,
        bstop,
        breset,
        constant_ic_items,
        constant_ic_textboxes,
        equation_items,
    )
end