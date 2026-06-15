# src/PerturbationPanel.jl

# ============================================================
# Perturbation panel
# ============================================================
#
# This file builds the local perturbation controls displayed
# under every plot.
#
# Each variable receives:
#
#     - center slider in [0, 1],
#     - Perturb button,
#     - constant/random switch,
#     - width numeric field with +/- buttons,
#     - height numeric field with +/- buttons.
#
# Moving any control stops the simulation and updates a gray
# preview curve.
#
# ============================================================


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
        # Assumes normalized domain [0, 1].
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
    # Small custom numeric field:
    #
    #     label
    #     [-] [value] [+]

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


function clear_preview_observable!(
    preview_obs::Observable{Vector{Float64}},
    N::Int,
)
    preview_obs[] = fill(NaN, N)
    return nothing
end


function update_local_perturbation_preview!(
    app::AppState,
    preview_obs::Observable{Vector{Float64}},
    variable::Int,
    center::Float64,
    width::Union{Nothing, Float64},
    height::Union{Nothing, Float64},
    random_mode::Bool,
)
    # Stop simulation and draw gray preview of the perturbation.

    stop_worker!(app; wait = true)

    lock(app.simlock)

    try
        if variable > app.sim.model.nvars
            return nothing
        end

        if width === nothing || height === nothing || width <= 0
            clear_preview_observable!(preview_obs, app.sim.N)
            return nothing
        end

        y = copy(app.sim.integrator_ref[].u)
        U = reshape(y, app.sim.N, app.sim.model.nvars)

        actual = copy(U[:, variable])
        preview = fill(NaN, app.sim.N)

        mask = local_perturbation_mask(
            app.sim.x,
            center,
            width;
            boundary_condition = app.sim.boundary_condition,
        )

        preview_height = random_mode ? 0.5 * height : height

        preview[mask] .= actual[mask] .+ preview_height

        preview_obs[] = preview

        rescale_axis_from_actual_and_preview!(
            app.plot_panel.axes[variable],
            actual,
            preview,
        )

    finally
        unlock(app.simlock)
    end

    return nothing
end


function build_perturbation_controls!(
    grid::GridLayout,
    app::AppState,
    preview_obs::Observable{Vector{Float64}};
    variable::Int,
)
    # Controls under one plot.

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

    width_grid = GridLayout()
    height_grid = GridLayout()

    grid[1:2, 4] = width_grid
    grid[1:2, 5] = height_grid

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

    append!(
        ui_items,
        Any[
            center_label,
            center_slider,
            perturb_button,
            mode_label,
            random_toggle,
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
        center = Float64(center_slider.value[])
        width = textbox_float_value(width_spinner.textbox; default = nothing)
        height = textbox_float_value(height_spinner.textbox; default = nothing)
        random_mode = random_toggle.active[]

        update_local_perturbation_preview!(
            app,
            preview_obs,
            variable,
            center,
            width,
            height,
            random_mode,
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

    on(perturb_button.clicks) do _
        center = Float64(center_slider.value[])
        width = textbox_float_value(width_spinner.textbox; default = nothing)
        height = textbox_float_value(height_spinner.textbox; default = nothing)
        random_mode = random_toggle.active[]

        if width === nothing || height === nothing || width <= 0
            @warn "Invalid perturbation parameters."
            return nothing
        end

        apply_local_perturbation_app!(
            app;
            variable = variable,
            center = center,
            width = width,
            height = height,
            random_mode = random_mode,
        )
    end

    return ui_items
end