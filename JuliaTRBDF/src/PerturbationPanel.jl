# src/PerturbationPanel.jl

# ============================================================
# Perturbation panel
# ============================================================


mutable struct PerturbationControlState
    variable::Int
    center_slider::Any
    random_toggle::Any
    set_toggle::Any
    width_textbox::Any
    height_textbox::Any
    preview_obs::Observable{Vector{Float64}}
    increment::Vector{Float64}
    has_valid_preview::Bool
end


function local_perturbation_mask(
    x::AbstractVector{<:Real},
    center::Real,
    width::Real;
    boundary_condition::Symbol,
)
    width > 0 ||
        return falses(length(x))

    half_width = width / 2

    if boundary_condition == :periodic
        return [min(abs(xi - center), 1.0 - abs(xi - center)) <= half_width for xi in x]
    else
        return [abs(xi - center) <= half_width for xi in x]
    end
end


function textbox_float_value(textbox; default = nothing)
    value = tryparse(Float64, textbox.stored_string[])

    if value === nothing
        return default
    end

    return value
end


function set_textbox_float!(textbox, value::Real; digits::Int = 4)
    textbox.stored_string[] = string(round(Float64(value); digits = digits))
    return nothing
end


function build_float_spinner!(
    grid::GridLayout;
    label::String,
    default::Float64,
    step::Float64,
    min_value::Float64 = -Inf,
    max_value::Float64 = Inf,
)
    label_block = Label(
        grid[1, 1:3],
        label,
        tellwidth = false,
    )

    minus_button = Button(
        grid[2, 1],
        label = "-",
        tellwidth = false,
    )

    textbox = Textbox(
        grid[2, 2],
        stored_string = string(default),
        tellwidth = false,
    )

    plus_button = Button(
        grid[2, 3],
        label = "+",
        tellwidth = false,
    )

    on(minus_button.clicks) do _
        old_value = textbox_float_value(textbox; default = default)
        new_value = clamp(old_value - step, min_value, max_value)
        set_textbox_float!(textbox, new_value)
    end

    on(plus_button.clicks) do _
        old_value = textbox_float_value(textbox; default = default)
        new_value = clamp(old_value + step, min_value, max_value)
        set_textbox_float!(textbox, new_value)
    end

    return (;
        label_block,
        minus_button,
        textbox,
        plus_button,
    )
end


function update_local_perturbation_preview!(
    app::AppState,
    state::PerturbationControlState;
    stop_simulation::Bool = true,
)
    # Update gray dashed preview curve.
    #
    # The preview is a full curve on the whole interval:
    #
    #     preview = current solution + perturbation increment.
    #
    # Outside the perturbation support, increment = 0, so the gray line
    # coincides with the original curve.

    if stop_simulation
        stop_worker!(app; wait = true)
    end

    success = false

    lock(app.simlock)

    try
        variable = state.variable

        if variable > app.sim.model.nvars
            state.has_valid_preview = false
            return false
        end

        width = textbox_float_value(state.width_textbox; default = nothing)
        height = textbox_float_value(state.height_textbox; default = nothing)

        if width === nothing || height === nothing || width <= 0
            N = app.sim.N
            state.preview_obs[] = fill(NaN, N)
            state.increment = zeros(Float64, N)
            state.has_valid_preview = false
            return false
        end

        center = Float64(state.center_slider.value[])
        random_mode = state.random_toggle.active[]
        set_mode = state.set_toggle.active[]

        y = copy(app.sim.integrator_ref[].u)
        U = reshape(y, app.sim.N, app.sim.model.nvars)

        actual = copy(U[:, variable])

        mask = local_perturbation_mask(
            app.sim.x,
            center,
            width;
            boundary_condition = app.sim.boundary_condition,
        )

        increment = zeros(Float64, app.sim.N)

        for i in eachindex(mask)
            if mask[i]
                value = random_mode ? height * rand() : height

                if set_mode
                    # Absolute mode:
                    #
                    #     preview[i] = value
                    #
                    # Since the runtime applies increments, we store
                    #
                    #     increment[i] = value - actual[i].
                    increment[i] = value - actual[i]
                else
                    # Perturbation mode:
                    #
                    #     preview[i] = actual[i] + value
                    increment[i] = value
                end
            end
        end

        preview = actual .+ increment

        state.increment = increment
        state.preview_obs[] = preview
        state.has_valid_preview = true

        rescale_axis_from_actual_and_preview!(
            app.plot_panel.axes[variable],
            actual,
            preview,
        )

        success = true

    finally
        unlock(app.simlock)
    end

    return success
end


function update_all_perturbation_previews!(
    app::AppState;
    stop_simulation::Bool = false,
)
    # Show gray dashed preview lines for all variables.
    #
    # This should be called after the simulation is stopped.

    if stop_simulation
        stop_worker!(app; wait = true)
    end

    for state in app.plot_panel.perturbation_controls
        update_local_perturbation_preview!(
            app,
            state;
            stop_simulation = false,
        )
    end

    return nothing
end


function build_perturbation_controls!(
    grid::GridLayout,
    app::AppState,
    preview_obs::Observable{Vector{Float64}};
    variable::Int,
)
    ui_items = Any[]

    center_label = Label(
        grid[1, 1],
        "center",
        tellwidth = false,
    )

    center_slider = Slider(
        grid[2, 1],
        range = 0.0:0.01:1.0,
        startvalue = 0.5,
        tellwidth = false,
    )

    perturb_button = Button(
        grid[2, 2],
        label = "Perturb",
        tellwidth = false,
    )

    mode_label = Label(
        grid[1, 3],
        "constant",
        tellwidth = false,
    )

    random_toggle = Toggle(
        grid[2, 3],
        active = false,
        tellwidth = false,
    )

    set_mode_label = Label(
        grid[1, 4],
        "perturb",
        tellwidth = false,
    )

    set_toggle = Toggle(
        grid[2, 4],
        active = false,
        tellwidth = false,
    )

    width_grid = GridLayout()
    height_grid = GridLayout()

    grid[1:2, 5] = width_grid
    grid[1:2, 6] = height_grid

    width_spinner = build_float_spinner!(
        width_grid;
        label = "width",
        default = 0.05,
        step = 0.01,
        min_value = 1e-6,
    )

    height_spinner = build_float_spinner!(
        height_grid;
        label = "height",
        default = 1.0,
        step = 0.1,
    )

    state = PerturbationControlState(
        variable,
        center_slider,
        random_toggle,
        set_toggle,
        width_spinner.textbox,
        height_spinner.textbox,
        preview_obs,
        zeros(Float64, app.sim.N),
        false,
    )

    append!(
        ui_items,
        Any[
            center_label,
            center_slider,
            perturb_button,
            mode_label,
            random_toggle,
            set_mode_label,
            set_toggle,
            width_grid,
            height_grid,
            width_spinner.label_block,
            width_spinner.minus_button,
            width_spinner.textbox,
            width_spinner.plus_button,
            height_spinner.label_block,
            height_spinner.minus_button,
            height_spinner.textbox,
            height_spinner.plus_button,
        ],
    )

    function update_preview()
        update_local_perturbation_preview!(
            app,
            state;
            stop_simulation = true,
        )

        return nothing
    end

    on(center_slider.value) do _
        update_preview()
    end

    on(width_spinner.textbox.stored_string) do _
        update_preview()
    end

    on(height_spinner.textbox.stored_string) do _
        update_preview()
    end

    on(random_toggle.active) do is_random
        mode_label.text[] = is_random ? "random" : "constant"
        update_preview()
    end

    on(set_toggle.active) do is_set
        set_mode_label.text[] = is_set ? "set" : "perturb"
    update_preview()
    end

    on(perturb_button.clicks) do _
        if !state.has_valid_preview
            ok = update_local_perturbation_preview!(
                app,
                state;
                stop_simulation = true,
            )

            ok || return nothing
        end

        apply_local_perturbation_increment_app!(
            app;
            variable = state.variable,
            increment = copy(state.increment),
        )
    end

    return (;
        ui_items,
        state,
    )
end