# src/PerturbationPanel.jl

# ============================================================
# Mouse-driven perturbation controls
# ============================================================


mutable struct PerturbationControlState
    random_mode::Observable{Bool}
    absolute_mode::Observable{Bool}
    width_textbox::Any
    height_textbox::Any
    width_value::Float64
    updating_width_textbox::Bool
    relative_axis_scaled::Bool
    active_segment::Int
    active_variable::Int
    center::Float64
    mouse_height::Float64
    increment::Vector{Float64}
    has_valid_preview::Bool
end


function active_perturbation_state(app::AppState)
    isempty(app.plot_panel.perturbation_controls) && return nothing

    state = first(app.plot_panel.perturbation_controls)
    return state isa PerturbationControlState ? state : nothing
end


function set_perturbation_random_mode!(app::AppState, random_mode::Bool)
    state = active_perturbation_state(app)
    state === nothing && return false
    state.random_mode[] = random_mode

    return true
end


function set_perturbation_absolute_mode!(app::AppState, absolute_mode::Bool)
    state = active_perturbation_state(app)
    state === nothing && return false
    state.absolute_mode[] = absolute_mode

    return true
end


function toggle_perturbation_random_mode!(app::AppState)
    state = active_perturbation_state(app)
    state === nothing && return false

    return set_perturbation_random_mode!(app, !state.random_mode[])
end


function toggle_perturbation_absolute_mode!(app::AppState)
    state = active_perturbation_state(app)
    state === nothing && return false

    return set_perturbation_absolute_mode!(app, !state.absolute_mode[])
end


function set_perturbation_width_value!(
    state::PerturbationControlState,
    width::Real,
)
    state.width_value = Float64(width)
    formatted_width = @sprintf("%.2f", state.width_value)

    state.updating_width_textbox = true

    try
        state.width_textbox.displayed_string[] = formatted_width
        state.width_textbox.stored_string[] = formatted_width
    finally
        state.updating_width_textbox = false
    end

    return nothing
end


function set_perturbation_width_value!(app::AppState, width::Real)
    state = active_perturbation_state(app)
    state === nothing && return false
    set_perturbation_width_value!(state, width)

    return true
end


function set_perturbation_height_value!(
    state::PerturbationControlState,
    height::Real,
)
    formatted_height = string(Float64(height))
    state.height_textbox.displayed_string[] = formatted_height
    state.height_textbox.stored_string[] = formatted_height

    return nothing
end


function set_perturbation_height_value!(app::AppState, height::Real)
    state = active_perturbation_state(app)
    state === nothing && return false
    set_perturbation_height_value!(state, height)

    return true
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
        period = (last(x) - first(x)) + (x[2] - x[1])

        return [min(abs(xi - center), period - abs(xi - center)) <= half_width for xi in x]
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


function perturbation_control_values(app::AppState)
    state = active_perturbation_state(app)
    state === nothing && return (width = 0.05, height = 0.0)

    height = textbox_float_value(state.height_textbox; default = 0.0)
    return (
        width = state.width_value,
        height = height === nothing ? 0.0 : height,
    )
end


function simulation_is_stopped(app::AppState)
    return !app.worker_running[] &&
           !app.running[] &&
           !app.synchronization_running[]
end


function simulation_domain_position(
    app::AppState,
    segment::Int,
    displayed_position::Real,
)
    sim = app.simulations[segment]
    displayed_length =
        segment_base_length(app, segment) *
        app.plot_panel.domain_length_scale
    fraction = displayed_length > 0 ?
        clamp(Float64(displayed_position) / displayed_length, 0.0, 1.0) :
        0.0
    xmin = first(sim.x)

    return xmin + fraction * simulation_domain_length(sim)
end


function simulation_domain_width(
    relative_width::Real,
    sim::SimulationState,
)
    return Float64(relative_width) * simulation_domain_length(sim)
end


function nearest_grid_index(
    x::AbstractVector{<:Real},
    value::Real,
)
    index = searchsortedfirst(x, value)

    if index <= 1
        return 1
    elseif index > length(x)
        return length(x)
    end

    left_index = index - 1

    if abs(x[index] - value) < abs(value - x[left_index])
        return index
    else
        return left_index
    end
end


function make_mouse_perturbation_increment(
    actual::AbstractVector{<:Real},
    mask::AbstractVector{Bool},
    x::AbstractVector{<:Real},
    center::Float64,
    mouse_height::Float64,
    absolute_height::Float64,
    random_mode::Bool,
    absolute_mode::Bool,
)
    increment = zeros(Float64, length(actual))

    if absolute_mode
        for i in eachindex(mask)
            if mask[i]
                target = random_mode ? absolute_height * rand() : absolute_height
                increment[i] = target - actual[i]
            end
        end
    else
        center_index = nearest_grid_index(x, center)
        relative_height = mouse_height - actual[center_index]

        for i in eachindex(mask)
            if mask[i]
                increment[i] = random_mode ? relative_height * rand() : relative_height
            end
        end
    end

    return increment
end


function clear_other_perturbation_previews!(
    panel::PlotPanel,
    active_segment::Int,
    active_variable::Int,
)
    for segment in eachindex(panel.segment_preview_observables)
        for variable in eachindex(panel.segment_preview_observables[segment])
            if segment != active_segment || variable != active_variable
                preview = panel.segment_preview_observables[segment][variable]
                N = length(preview[])
                preview[] = fill(NaN, N)
            end
        end
    end

    return nothing
end


function update_mouse_perturbation_preview!(
    app::AppState,
    state::PerturbationControlState;
    segment::Int,
    variable::Int,
    center::Float64,
    mouse_height::Float64,
)
    if !simulation_is_stopped(app)
        clear_perturbation_previews!(app.plot_panel)
        return false
    end

    width = state.width_value

    if !isfinite(width) || width <= 0 || width > 1
        clear_perturbation_previews!(app.plot_panel)
        return false
    end

    absolute_height = 0.0

    if state.absolute_mode[]
        parsed_height = textbox_float_value(
            state.height_textbox;
            default = nothing,
        )

        if parsed_height === nothing || !isfinite(parsed_height)
            clear_perturbation_previews!(app.plot_panel)
            return false
        end

        absolute_height = parsed_height
    end

    success = false

    lock(app.simlock)

    try
        if !simulation_is_stopped(app) ||
           segment < 1 ||
           segment > length(app.simulations) ||
           variable < 1 ||
           variable > app.sim.model.nvars
            clear_perturbation_previews!(app.plot_panel)
            return false
        end

        sim = app.simulations[segment]
        y = copy(sim.integrator_ref[].u)
        U = reshape(y, sim.N, sim.model.nvars)
        actual = copy(U[:, variable])

        simulation_center = simulation_domain_position(
            app,
            segment,
            center,
        )

        simulation_width = simulation_domain_width(
            width,
            sim,
        )

        mask = local_perturbation_mask(
            sim.x,
            simulation_center,
            simulation_width;
            boundary_condition = sim.boundary_condition,
        )

        if !any(mask)
            mask[nearest_grid_index(sim.x, simulation_center)] = true
        end

        increment = make_mouse_perturbation_increment(
            actual,
            mask,
            sim.x,
            simulation_center,
            mouse_height,
            absolute_height,
            state.random_mode[],
            state.absolute_mode[],
        )

        preview = actual .+ increment

        clear_other_perturbation_previews!(
            app.plot_panel,
            segment,
            variable,
        )

        state.active_segment = segment
        state.active_variable = variable
        state.center = center
        state.mouse_height = mouse_height
        state.increment = increment
        state.has_valid_preview = true

        app.plot_panel.segment_preview_observables[segment][variable][] = preview

        if state.absolute_mode[]
            rescale_solution_variable_axes!(app.plot_panel, variable)
        end

        success = true

    finally
        unlock(app.simlock)
    end

    return success
end


function change_perturbation_width_from_scroll!(
    app::AppState,
    state::PerturbationControlState,
    segment::Int,
    scroll_delta::Real,
)
    delta = clamp(Float64(scroll_delta), -10.0, 10.0)

    iszero(delta) &&
        return false

    current_width = state.width_value
    sim = app.simulations[segment]
    simulation_length = simulation_domain_length(sim)

    if !isfinite(current_width) || current_width <= 0
        current_width = 0.05
    end

    minimum_width = min(
        1.0,
        max(sim.dx / simulation_length, eps(1.0)),
    )

    new_width = clamp(
        current_width * 1.1^delta,
        minimum_width,
        1.0,
    )

    set_perturbation_width_value!(state, new_width)

    return true
end


function change_perturbation_height_from_scroll!(
    state::PerturbationControlState,
    scroll_delta::Real,
)
    state.absolute_mode[] ||
        return false

    direction = sign(Float64(scroll_delta))

    iszero(direction) &&
        return false

    current_height = textbox_float_value(
        state.height_textbox;
        default = 0.0,
    )

    if !isfinite(current_height)
        current_height = 0.0
    end

    new_height = current_height + direction
    set_perturbation_height_value!(state, new_height)

    return true
end


function amplify_relative_perturbation_from_scroll!(
    app::AppState,
    state::PerturbationControlState,
    scroll_delta::Real,
)
    state.absolute_mode[] &&
        return false

    simulation_is_stopped(app) ||
        return false

    state.has_valid_preview ||
        return false

    direction = sign(Float64(scroll_delta))

    iszero(direction) &&
        return false

    segment = state.active_segment
    variable = state.active_variable

    if segment < 1 || segment > length(app.simulations) ||
       variable < 1 || variable > app.sim.model.nvars
        return false
    end

    scale_factor = 1.5^direction

    lock(app.simlock)

    try
        simulation_is_stopped(app) ||
            return false

        sim = app.simulations[segment]
        y = copy(sim.integrator_ref[].u)
        U = reshape(y, sim.N, sim.model.nvars)
        actual = copy(U[:, variable])

        increment = state.increment .* scale_factor
        preview = actual .+ increment

        state.increment = increment
        state.relative_axis_scaled = true
        app.plot_panel.segment_preview_observables[segment][variable][] = preview

        rescale_solution_variable_axes!(app.plot_panel, variable)
    finally
        unlock(app.simlock)
    end

    return true
end


function restore_relative_axis_from_solution!(
    app::AppState,
    state::PerturbationControlState,
    ;
    force::Bool = false,
)
    (force || state.relative_axis_scaled) ||
        return false

    rescale_solution_axes!(
        app.plot_panel,
        app.simulations;
        include_previews = false,
    )

    state.relative_axis_scaled = false

    return true
end


function register_perturbation_scroll_handlers!(
    app::AppState,
    segment_axes::Vector{Vector{Axis}},
    state::PerturbationControlState,
)
    for (segment, axes) in enumerate(segment_axes)
        for axis in axes
            on(events(axis.scene).scroll, priority = 20) do scroll
            if !is_mouseinside(axis.scene)
                return Consume(false)
            end

            scroll_delta = iszero(scroll[2]) ? scroll[1] : scroll[2]
            control_pressed = ispressed(
                axis.scene,
                Keyboard.left_control | Keyboard.right_control,
            )

            if control_pressed
                if state.absolute_mode[]
                    change_perturbation_height_from_scroll!(
                        state,
                        scroll_delta,
                    )
                else
                    amplify_relative_perturbation_from_scroll!(
                        app,
                        state,
                        scroll_delta,
                    )
                end
            else
                change_perturbation_width_from_scroll!(
                    app,
                    state,
                    segment,
                    scroll_delta,
                )
            end

            return Consume(true)
            end
        end
    end

    if !isempty(segment_axes) && !isempty(first(segment_axes))
        keyboard_scene = first(first(segment_axes)).scene

        on(events(keyboard_scene).keyboardbutton, priority = 20) do event
            is_control_release =
                event.action == Keyboard.release &&
                event.key in (
                    Keyboard.left_control,
                    Keyboard.right_control,
                )

            control_still_pressed = ispressed(
                keyboard_scene,
                Keyboard.left_control | Keyboard.right_control,
            )

            if is_control_release &&
               !control_still_pressed &&
               !state.absolute_mode[] &&
               restore_relative_axis_from_solution!(app, state)
                if !state.absolute_mode[]
                    update_perturbation_preview_from_mouse!(
                        app,
                        state,
                    )
                end
            end

            return Consume(false)
        end
    end

    return nothing
end


function update_perturbation_preview_from_mouse!(
    app::AppState,
    state::PerturbationControlState,
)
    if !simulation_is_stopped(app)
        clear_perturbation_previews!(app.plot_panel)
        return false
    end

    for segment in eachindex(app.plot_panel.segment_axes)
        for variable in eachindex(app.plot_panel.segment_axes[segment])
            axis = app.plot_panel.segment_axes[segment][variable]

            if is_mouseinside(axis.scene)
                position = mouseposition(axis.scene)

                return update_mouse_perturbation_preview!(
                    app,
                    state;
                    segment = segment,
                    variable = variable,
                    center = Float64(position[1]),
                    mouse_height = Float64(position[2]),
                )
            end
        end
    end

    clear_perturbation_previews!(app.plot_panel)

    return false
end


function update_all_perturbation_previews!(
    app::AppState;
    stop_simulation::Bool = false,
)
    if stop_simulation
        stop_worker!(app; wait = true)
    end

    isempty(app.plot_panel.perturbation_controls) &&
        return nothing

    state = first(app.plot_panel.perturbation_controls)

    update_perturbation_preview_from_mouse!(
        app,
        state,
    )

    return nothing
end


function register_mouse_perturbation_handlers!(
    app::AppState,
    segment_axes::Vector{Vector{Axis}},
    state::PerturbationControlState,
)
    for (segment, axes) in enumerate(segment_axes)
        for (variable, axis) in enumerate(axes)
        on(events(axis.scene).mouseposition, priority = 10) do _
            if is_mouseinside(axis.scene)
                position = mouseposition(axis.scene)

                update_mouse_perturbation_preview!(
                    app,
                    state;
                    segment = segment,
                    variable = variable,
                    center = Float64(position[1]),
                    mouse_height = Float64(position[2]),
                )

            elseif state.active_segment == segment &&
                   state.active_variable == variable
                clear_perturbation_previews!(app.plot_panel)
            end

            return Consume(false)
        end

        on(events(axis.scene).mousebutton, priority = 10) do event
            is_left_press =
                event.button == Mouse.left &&
                event.action == Mouse.press

            if !is_left_press || !is_mouseinside(axis.scene)
                return Consume(false)
            end

            if !simulation_is_stopped(app)
                stop_worker!(app; wait = true)

                position = mouseposition(axis.scene)

                update_mouse_perturbation_preview!(
                    app,
                    state;
                    segment = segment,
                    variable = variable,
                    center = Float64(position[1]),
                    mouse_height = Float64(position[2]),
                )

                return Consume(true)
            end

            if state.active_segment != segment ||
               state.active_variable != variable ||
               !state.has_valid_preview
                position = mouseposition(axis.scene)

                update_mouse_perturbation_preview!(
                    app,
                    state;
                    segment = segment,
                    variable = variable,
                    center = Float64(position[1]),
                    mouse_height = Float64(position[2]),
                )
            end

            if state.active_segment == segment &&
               state.active_variable == variable &&
               state.has_valid_preview
                apply_local_perturbation_increment_app!(
                    app;
                    segment = segment,
                    variable = variable,
                    increment = copy(state.increment),
                )

                return Consume(true)
            end

            return Consume(false)
            end
        end
    end

    return nothing
end


function build_perturbation_controls!(
    grid::GridLayout,
    app::AppState;
    segment_axes::Vector{Vector{Axis}},
)
    random_mode = Observable(true)
    absolute_mode = Observable(false)
    initial_width = 0.05

    color_inactive = :lightgray
    color_active = :skyblue

    random_button = Button(
        grid[1, 1],
        label = lift(random_mode) do is_random
            is_random ? "Random" : "Constant"
        end,
        buttoncolor = lift(random_mode) do is_random
            is_random ? color_active : color_inactive
        end,
        tellwidth = false,
    )

    mode_button = Button(
        grid[1, 2],
        label = lift(absolute_mode) do is_absolute
            is_absolute ? "Absolute" : "Relative"
        end,
        buttoncolor = lift(absolute_mode) do is_absolute
            is_absolute ? color_active : color_inactive
        end,
        tellwidth = false,
    )

    width_label = Label(
        grid[1, 3],
        "Width",
        tellwidth = false,
    )

    width_textbox = Textbox(
        grid[1, 4],
        stored_string = @sprintf("%.2f", initial_width),
        width = 70,
        tellwidth = false,
    )

    height_label = Label(
        grid[1, 5],
        "Height",
        tellwidth = false,
    )

    height_textbox = Textbox(
        grid[1, 6],
        stored_string = "0.0",
        width = 70,
        tellwidth = false,
    )

    state = PerturbationControlState(
        random_mode,
        absolute_mode,
        width_textbox,
        height_textbox,
        initial_width,
        false,
        false,
        0,
        0,
        NaN,
        NaN,
        zeros(Float64, app.sim.N),
        false,
    )

    function update_height_visibility!()
        is_visible = app.show_embedded_perturbation_controls && absolute_mode[]

        for block in (height_label, height_textbox)
            block.blockscene.visible[] = is_visible

            if hasproperty(block, :scene) && isdefined(block, :scene)
                block.scene.visible[] = is_visible
            end
        end

        return nothing
    end

    function update_preview_from_current_mouse!()
        update_perturbation_preview_from_mouse!(
            app,
            state,
        )

        return nothing
    end

    on(random_button.clicks) do _
        random_mode[] = !random_mode[]
        return nothing
    end

    on(random_mode) do _
        update_preview_from_current_mouse!()
        return nothing
    end

    on(mode_button.clicks) do _
        absolute_mode[] = !absolute_mode[]
        return nothing
    end

    on(absolute_mode) do _
        restore_relative_axis_from_solution!(
            app,
            state;
            force = !absolute_mode[],
        )

        update_height_visibility!()
        update_preview_from_current_mouse!()
        return nothing
    end

    on(width_textbox.stored_string) do value
        if !state.updating_width_textbox
            parsed_width = tryparse(Float64, value)
            state.width_value = parsed_width === nothing ? NaN : parsed_width
        end

        update_preview_from_current_mouse!()
        return nothing
    end

    on(height_textbox.stored_string) do _
        absolute_mode[] && update_preview_from_current_mouse!()
        return nothing
    end

    update_height_visibility!()

    register_mouse_perturbation_handlers!(
        app,
        segment_axes,
        state,
    )

    ui_items = Any[
        random_button,
        mode_button,
        width_label,
        width_textbox,
        height_label,
        height_textbox,
    ]

    colgap!(grid, 10)
    colsize!(grid, 1, Fixed(110))
    colsize!(grid, 2, Fixed(110))
    colsize!(grid, 3, Fixed(55))
    colsize!(grid, 4, Fixed(80))
    colsize!(grid, 5, Fixed(55))
    colsize!(grid, 6, Fixed(80))

    return (;
        ui_items,
        state,
    )
end
