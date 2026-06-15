# src/PerturbationPanel.jl

# ============================================================
# Perturbation panel
# ============================================================


mutable struct PerturbationControlState
    variable::Int
    center_slider::Any
    random_mode::Observable{Bool}
    set_mode::Observable{Bool}
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


function update_local_perturbation_preview!(
    app::AppState,
    state::PerturbationControlState;
    stop_simulation::Bool = true,
)
    # Update gray dashed preview curve.
    #
    # The preview is a full curve:
    #
    #     preview = current solution + perturbation increment.
    #
    # Outside the perturbation support, increment = 0, so the gray dashed
    # line coincides with the current solution.

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
        random_mode = state.random_mode[]
        set_mode = state.set_mode[]

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
                    # Since the runtime applies increments, store
                    #
                    #     increment[i] = value - actual[i].
                    increment[i] = value - actual[i]
                else
                    # Perturbation mode:
                    #
                    #     preview[i] = actual[i] + value.
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

    random_mode = Observable(false)
    set_mode = Observable(false)

    # Definujemy kolory dla stanów guzików dwustanowych
    color_inactive = :lightgray
    color_active = :skyblue    # Ładny, wyróżniający się kolor dla stanu "włączony"

    # 1. Główny Slider na całą wolną przestrzeń
    center_slider = Slider(
        grid[1, 1],
        range = 0.0:0.01:1.0,
        startvalue = 0.5,
    )

    # 2. Podsiatka na skurczone kontrolki
    ctrl_grid = grid[1, 2] = GridLayout()

    # Kolumna 1: Przycisk akcji
    perturb_button = Button(
        ctrl_grid[1, 1],
        label = "Perturb",
    )

    # Kolumna 2: Guzik Constant / Random
    random_button = Button(
        ctrl_grid[1, 2],
        label = "Constant",
        buttoncolor = color_inactive,
    )

    # Kolumna 3: Guzik Relative / Absolute
    set_button = Button(
        ctrl_grid[1, 3],
        label = "Relative",
        buttoncolor = color_inactive,
    )

    # Kolumna 4 i 5: Pole Width
    width_label = Label(
        ctrl_grid[1, 4],
        "Width",
    )
    width_textbox = Textbox(
        ctrl_grid[1, 5],
        stored_string = "0.05",
        width = 50,
    )

    # Kolumna 6 i 7: Pole Height
    height_label = Label(
        ctrl_grid[1, 6],
        "Height",
    )
    height_textbox = Textbox(
        ctrl_grid[1, 7],
        stored_string = "1.0",
        width = 50,
    )

    # 3. Precyzyjne zarządzanie odstępami (Gaps)
    colgap!(ctrl_grid, 4)          # Standardowy, bardzo ciasny odstęp między elementami (4px)
    
    colgap!(ctrl_grid, 1, 20)      # Odstęp PO "Perturb" a PRZED "Constant/Random"
    colgap!(ctrl_grid, 3, 20)      # Odstęp PO "Relative/Absolute" a PRZED "Width"
    colgap!(ctrl_grid, 5, 20)      # Odstęp PO polu tekstowym width a PRZED "Height"
    
    colgap!(grid, 15)              # Odstęp między głównym sliderem a panelem kontrolek

    # Ograniczenia wielkości głównych kolumn
    colsize!(grid, 1, Auto(true))
    colsize!(grid, 2, Auto(false))

    state = PerturbationControlState(
        variable,
        center_slider,
        random_mode,
        set_mode,
        width_textbox,
        height_textbox,
        preview_obs,
        zeros(Float64, app.sim.N),
        false,
    )

    append!(
        ui_items,
        Any[
            center_slider,
            perturb_button,
            random_button,
            set_button,
            width_label,
            width_textbox,
            height_label,
            height_textbox,
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

    on(width_textbox.stored_string) do _
        update_preview()
    end

    on(height_textbox.stored_string) do _
        update_preview()
    end

    # Dynamiczny guzik: Constant / Random
    on(random_button.clicks) do _
        random_mode[] = !random_mode[]
        
        # Zmiana tekstu i koloru tła
        random_button.label[] = random_mode[] ? "Random" : "Constant"
        random_button.buttoncolor = random_mode[] ? color_active : color_inactive
        
        update_preview()
    end

    # Dynamiczny guzik: Relative / Absolute
    on(set_button.clicks) do _
        set_mode[] = !set_mode[]
        
        # Zmiana tekstu i koloru tła
        set_button.label[] = set_mode[] ? "Absolute" : "Relative"
        set_button.buttoncolor = set_mode[] ? color_active : color_inactive
        
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