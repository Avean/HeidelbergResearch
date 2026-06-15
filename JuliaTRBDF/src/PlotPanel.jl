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
#     - one local perturbation control row below every plot,
#     - gray preview curves for perturbations.
#
# ============================================================


function empty_plot_panel()
    return PlotPanel(
        Axis[],
        Observable{Vector{Float64}}[],
        Observable{Vector{Float64}}[],
        Any[],
    )
end


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

    if ymin == ymax
        pad = max(1.0, abs(ymin)) * 0.1
    else
        pad = 0.05 * (ymax - ymin)
    end

    ylims!(ax, ymin - pad, ymax + pad)

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


function build_plot_panel!(
    grid::GridLayout,
    app::AppState;
    title_obs = nothing,
)
    sim = app.sim
    model = sim.model
    U = solution_matrix(sim)

    axes = Axis[]
    observables = Observable{Vector{Float64}}[]
    preview_observables = Observable{Vector{Float64}}[]
    ui_items = Any[]

    for j in 1:model.nvars
        plot_row = 2j - 1
        control_row = 2j

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

        y_obs = Observable(copy(U[:, j]))

        lines!(
            ax,
            sim.x,
            y_obs,
            linewidth = 2,
        )

        preview_obs = Observable(fill(NaN, sim.N))

        lines!(
            ax,
            sim.x,
            preview_obs;
            color = (:gray, 0.35),
            linewidth = 5,
        )

        perturbation_grid = GridLayout()
        grid[control_row, 1] = perturbation_grid

        items = build_perturbation_controls!(
            perturbation_grid,
            app,
            preview_obs;
            variable = j,
        )

        append!(ui_items, items)

        push!(axes, ax)
        push!(observables, y_obs)
        push!(preview_observables, preview_obs)

        rowsize!(grid, control_row, Fixed(70))
    end

    panel = PlotPanel(
        axes,
        observables,
        preview_observables,
        ui_items,
    )

    rescale_solution_axes!(panel, sim)

    return panel
end


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


function clear_perturbation_preview!(
    panel::PlotPanel,
    variable::Int,
)
    if variable >= 1 && variable <= length(panel.preview_observables)
        N = length(panel.preview_observables[variable][])
        panel.preview_observables[variable][] = fill(NaN, N)
    end

    return nothing
end


function clear_perturbation_previews!(panel::PlotPanel)
    for j in eachindex(panel.preview_observables)
        clear_perturbation_preview!(panel, j)
    end

    return nothing
end


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
    empty!(panel.observables)
    empty!(panel.preview_observables)
    empty!(panel.ui_items)

    return nothing
end