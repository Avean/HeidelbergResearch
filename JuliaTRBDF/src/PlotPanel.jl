# src/PlotPanel.jl

# ============================================================
# Plot panel
# ============================================================
#
# Left application area.
#
# Contains:
#
#     - one plot for every model variable,
#     - one shared perturbation control row below the solution plots,
#     - gray dashed preview curves for local perturbations,
#     - optional spatial profile plots,
#     - Previous/Next selector for active spatial profile sets.
#
# ============================================================


# ============================================================
# Empty panel
# ============================================================

function empty_plot_panel()
    return PlotPanel(
        Axis[],
        Observable(Float64[]),
        1.0,
        Observable{Vector{Float64}}[],
        Observable{Vector{Float64}}[],
        Any[],
        Any[],
    )
end


function displayed_domain_coordinates(
    sim::SimulationState,
    domain_length_scale::Float64,
)
    xmin = first(sim.x)

    return @. xmin + domain_length_scale * (sim.x - xmin)
end


function simulation_domain_length(sim::SimulationState)
    interval_count =
        sim.boundary_condition == :periodic ? sim.N : sim.N - 1

    return interval_count * sim.dx
end


function set_plot_domain_scale!(
    panel::PlotPanel,
    sim::SimulationState,
    diffusion_scale::Real,
)
    scale = Float64(diffusion_scale)

    isfinite(scale) && scale > 0 ||
        error("Domain scale must be positive and finite.")

    domain_length_scale = sqrt(scale)
    displayed_x = displayed_domain_coordinates(
        sim,
        domain_length_scale,
    )

    panel.domain_length_scale = domain_length_scale
    panel.x_observable[] = displayed_x

    xmin = first(sim.x)
    xmax = xmin + domain_length_scale * simulation_domain_length(sim)

    for axis in panel.axes
        xlims!(axis, xmin, xmax)
    end

    clear_perturbation_previews!(panel)

    return nothing
end


# ============================================================
# Axis scaling
# ============================================================

function finite_values(v::AbstractVector)
    return collect(filter(isfinite, v))
end


function set_axis_y_limits_from_values!(
    ax::Axis,
    values::AbstractVector,
)
    finite = finite_values(values)

    isempty(finite) &&
        return nothing

    ymin, ymax = extrema(finite)
    ymin = 0.0
    ymean = mean(finite)

    # Minimal allowed y-range.
    # This prevents degenerate or almost-degenerate axis limits.
    min_range = 1e-3

    yrange = ymax - ymin

    if yrange < min_range
        center = 0.5 * (ymin + ymax)
        ymin = center - 0.5 * min_range
        ymax = center + 0.5 * min_range
    else
        pad = 0.05 * yrange
        ymin -= pad
        ymax += pad
    end

    ylims!(ax, ymin, max(ymax, 2*ymean))

    return nothing
end


function rescale_axis_from_actual_and_preview!(
    ax::Axis,
    actual::AbstractVector,
    preview::AbstractVector,
)
    values = Float64[]

    append!(values, finite_values(actual))
    append!(values, finite_values(preview))

    set_axis_y_limits_from_values!(ax, values)

    return nothing
end


function rescale_solution_axes!(
    panel::PlotPanel,
    sim::SimulationState,
)
    U = solution_matrix(sim)

    for j in 1:sim.model.nvars
        actual = U[:, j]
        preview = panel.preview_observables[j][]

        rescale_axis_from_actual_and_preview!(
            panel.axes[j],
            actual,
            preview,
        )
    end

    return nothing
end


function rescale_solution_axes_from_snapshot!(
    panel::PlotPanel,
    snapshot::SimulationSnapshot,
)
    U = solution_matrix_from_snapshot(snapshot)

    for j in 1:snapshot.nvars
        actual = U[:, j]
        preview = panel.preview_observables[j][]

        rescale_axis_from_actual_and_preview!(
            panel.axes[j],
            actual,
            preview,
        )
    end

    return nothing
end


# ============================================================
# Spatial profiles
# ============================================================

function evaluate_spatial_profile(
    sim::SimulationState,
    profile_name::String,
    profile_fun::Function,
)
    raw = profile_fun(sim.x, sim.params)

    y = if raw isa Number
        fill(Float64(raw), sim.N)
    else
        Float64.(collect(raw))
    end

    length(y) == sim.N ||
        error("Spatial profile $(profile_name) has wrong length.")

    return y
end


function build_spatial_profile_panel!(
    grid::GridLayout,
    app::AppState,
    axes::Vector{Axis},
    x_observable::Observable{Vector{Float64}},
    ui_items::Vector{Any};
    start_row::Int,
)
    sim = app.sim
    profile_sets = sim.model.spatial_profile_sets

    isempty(profile_sets) &&
        return nothing

    max_profiles = maximum(length(profile_set) for (_, profile_set) in profile_sets)

    # --------------------------------------------------------
    # Profile plots
    # --------------------------------------------------------

    profile_axes = Axis[]
    profile_observables = Observable{Vector{Float64}}[]
    profile_name_observables = Observable{String}[]

    for k in 1:max_profiles
        row = start_row + k - 1

        profile_name_obs = Observable("")

        ax = Axis(
            grid[row, 1],
            xlabel = k == max_profiles ? "x" : "",
            ylabel = profile_name_obs,
            title = profile_name_obs,
        )

        deactivate_interaction!(ax, :scrollzoom)

        y_obs = Observable(fill(NaN, sim.N))

        lines!(
            ax,
            x_observable,
            y_obs,
            linewidth = 4,
            color = :red,
            linestyle = :dash,
        )

        push!(axes, ax)
        push!(profile_axes, ax)
        push!(profile_observables, y_obs)
        push!(profile_name_observables, profile_name_obs)

        rowsize!(grid, row, Fixed(120))
    end

    # --------------------------------------------------------
    # Navigation row under the profile plots
    # --------------------------------------------------------

    nav_row = start_row + max_profiles

    nav_grid = GridLayout(
        tellwidth = false,
        tellheight = true,
    )

    grid[nav_row, 1] = nav_grid

    previous_button = Button(
        nav_grid[1, 1],
        label = "Previous",
        tellwidth = false,
    )

    set_label_obs = Observable("spatial profile: ")

    set_label = Label(
        nav_grid[1, 2],
        set_label_obs,
        tellwidth = false,
        halign = :center,
    )

    next_button = Button(
        nav_grid[1, 3],
        label = "Next",
        tellwidth = false,
    )

    colsize!(nav_grid, 1, Fixed(80))
    colsize!(nav_grid, 2, Fixed(190))
    colsize!(nav_grid, 3, Fixed(80))
    colgap!(nav_grid, 5)

    try
        nav_grid.halign = :center
    catch
    end

    rowsize!(grid, nav_row, Fixed(34))

    push!(ui_items, nav_grid)
    push!(ui_items, previous_button)
    push!(ui_items, set_label)
    push!(ui_items, next_button)

    # --------------------------------------------------------
    # Current profile set
    # --------------------------------------------------------

    current_set_index = Ref(
        _active_spatial_profile_set_index(
            sim.params,
            profile_sets,
        ),
    )

    function update_profile_set!()
        current_sim = app.sim
        current_profile_sets = current_sim.model.spatial_profile_sets

        isempty(current_profile_sets) &&
            return nothing

        current_set_index[] = clamp(
            current_set_index[],
            1,
            length(current_profile_sets),
        )

        set_name, profiles = current_profile_sets[current_set_index[]]

        set_label_obs[] = "spatial profile: $(set_name)"

        for k in 1:max_profiles
            row = start_row + k - 1

            if k <= length(profiles)
                profile_name, profile_fun = profiles[k]

                y = evaluate_spatial_profile(
                    current_sim,
                    profile_name,
                    profile_fun,
                )

                profile_name_observables[k][] = profile_name
                profile_observables[k][] = y

                rowsize!(grid, row, Fixed(120))

                set_axis_y_limits_from_values!(
                    profile_axes[k],
                    y,
                )
            else
                profile_name_observables[k][] = ""
                profile_observables[k][] = fill(NaN, current_sim.N)

                rowsize!(grid, row, Fixed(0))
            end
        end

        return nothing
    end

    function switch_profile_set!(new_index::Int)
        current_set_index[] = new_index

        set_active_spatial_profile_set_app!(
            app,
            new_index,
        )

        update_profile_set!()

        return nothing
    end

    on(previous_button.clicks) do _
        nsets = length(app.sim.model.spatial_profile_sets)

        nsets == 0 &&
            return nothing

        new_index =
            current_set_index[] == 1 ? nsets : current_set_index[] - 1

        switch_profile_set!(new_index)
    end

    on(next_button.clicks) do _
        nsets = length(app.sim.model.spatial_profile_sets)

        nsets == 0 &&
            return nothing

        new_index =
            current_set_index[] == nsets ? 1 : current_set_index[] + 1

        switch_profile_set!(new_index)
    end

    update_profile_set!()

    return nothing
end

# ============================================================
# Build plot panel
# ============================================================

function build_plot_panel!(
    grid::GridLayout,
    app::AppState;
    title_obs = nothing,
)
    sim = app.sim
    model = sim.model
    U = solution_matrix(sim)
    x_observable = Observable(copy(sim.x))

    axes = Axis[]
    observables = Observable{Vector{Float64}}[]
    preview_observables = Observable{Vector{Float64}}[]
    perturbation_controls = Any[]
    ui_items = Any[]

    # Keep the plots stretched to the full width of the main window.
    # The shared perturbation toolbar has fixed-width controls and must not
    # determine the width of the plot column.
    colsize!(grid, 1, Relative(1.0))
    rowgap!(grid, 0)

    for j in 1:model.nvars
        plot_row = j

        axis_title = if j == 1 && title_obs !== nothing
            title_obs
        else
            model.varnames[j]
        end

        ax = Axis(
            grid[plot_row, 1],
            xlabel = j == model.nvars ? "x" : "",
            ylabel = model.varnames[j],
            title = axis_title,
        )

        deactivate_interaction!(ax, :scrollzoom)

        y_obs = Observable(copy(U[:, j]))

        lines!(
            ax,
            x_observable,
            y_obs,
            linewidth = 2,
        )

        preview_obs = Observable(fill(NaN, sim.N))

        lines!(
            ax,
            x_observable,
            preview_obs;
            color = (:gray, 0.45),
            linewidth = 2,
            linestyle = :dash,
        )

        push!(axes, ax)
        push!(observables, y_obs)
        push!(preview_observables, preview_obs)

        rowsize!(grid, plot_row, Auto(false, 1.0))
    end

    perturbation_row = model.nvars + 1
    perturbation_grid = GridLayout()
    grid[perturbation_row, 1] = perturbation_grid

    perturbation_panel = build_perturbation_controls!(
        perturbation_grid,
        app;
        axes = axes,
    )

    push!(ui_items, perturbation_grid)
    append!(ui_items, perturbation_panel.ui_items)
    push!(perturbation_controls, perturbation_panel.state)

    rowsize!(grid, perturbation_row, Fixed(42))

    profile_start_row = perturbation_row + 1

    build_spatial_profile_panel!(
        grid,
        app,
        axes,
        x_observable,
        ui_items;
        start_row = profile_start_row,
    )

    register_perturbation_scroll_handlers!(
        app,
        axes,
        perturbation_panel.state,
    )

    panel = PlotPanel(
        axes,
        x_observable,
        1.0,
        observables,
        preview_observables,
        perturbation_controls,
        ui_items,
    )

    set_plot_domain_scale!(panel, sim, 1.0)
    rescale_solution_axes!(panel, sim)

    return panel
end


# ============================================================
# Refresh
# ============================================================

function solution_matrix_from_snapshot(snapshot::SimulationSnapshot)
    return reshape(snapshot.y, snapshot.N, snapshot.nvars)
end


function refresh_plot_panel!(
    panel::PlotPanel,
    sim::SimulationState,
)
    U = solution_matrix(sim)

    length(panel.observables) == sim.model.nvars ||
        error("Plot panel does not match the number of model variables.")

    for j in 1:sim.model.nvars
        panel.observables[j][] = copy(U[:, j])
    end

    rescale_solution_axes!(panel, sim)

    return nothing
end


function refresh_plot_panel_from_snapshot!(
    panel::PlotPanel,
    snapshot::SimulationSnapshot,
)
    U = solution_matrix_from_snapshot(snapshot)

    length(panel.observables) == snapshot.nvars ||
        error("Plot panel does not match the number of snapshot variables.")

    for j in 1:snapshot.nvars
        panel.observables[j][] = copy(U[:, j])
    end

    rescale_solution_axes_from_snapshot!(panel, snapshot)

    return nothing
end


# ============================================================
# Perturbation previews
# ============================================================

function clear_perturbation_preview!(
    panel::PlotPanel,
    variable::Int,
)
    if variable >= 1 && variable <= length(panel.preview_observables)
        N = length(panel.preview_observables[variable][])
        panel.preview_observables[variable][] = fill(NaN, N)

        if !isempty(panel.perturbation_controls)
            state = first(panel.perturbation_controls)

            if state.active_variable != variable
                return nothing
            end

            state.active_variable = 0
            state.center = NaN
            state.mouse_height = NaN
            state.increment = zeros(Float64, N)
            state.has_valid_preview = false
        end
    end

    return nothing
end


function clear_perturbation_previews!(panel::PlotPanel)
    for preview_obs in panel.preview_observables
        N = length(preview_obs[])
        preview_obs[] = fill(NaN, N)
    end

    if !isempty(panel.perturbation_controls)
        state = first(panel.perturbation_controls)
        N = isempty(panel.preview_observables) ? 0 : length(first(panel.preview_observables)[])

        state.active_variable = 0
        state.center = NaN
        state.mouse_height = NaN
        state.increment = zeros(Float64, N)
        state.has_valid_preview = false
    end

    return nothing
end


# ============================================================
# Clearing / rebuilding
# ============================================================

function delete_plot_panel_item!(item)
    try
        delete!(item)
    catch
        try
            item.visible = false
        catch
        end
    end

    return nothing
end


function clear_plot_panel!(panel::PlotPanel)
    for ax in panel.axes
        delete_plot_panel_item!(ax)
    end

    for item in panel.ui_items
        delete_plot_panel_item!(item)
    end

    empty!(panel.axes)
    panel.x_observable[] = Float64[]
    panel.domain_length_scale = 1.0
    empty!(panel.observables)
    empty!(panel.preview_observables)
    empty!(panel.perturbation_controls)
    empty!(panel.ui_items)

    return nothing
end
