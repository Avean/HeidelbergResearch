# src/PlotPanel.jl

# ============================================================
# Plot panel
# ============================================================
#
# This file contains plotting utilities for displaying the
# current solution.
#
# For a model with nvars variables, we create nvars plots:
#
#     u_1(x)
#     u_2(x)
#     ...
#     u_n(x)
#
# Each variable gets its own vertical axis and its own automatic
# y-axis scaling.
#
# ============================================================


function empty_plot_panel()
    # Create an empty plot panel.
    #
    # The actual Makie axes and observables are created later by
    # build_plot_panel!.

    return PlotPanel(
        Axis[],
        Observable{Vector{Float64}}[],
    )
end


function build_plot_panel!(
    grid::GridLayout,
    sim::SimulationState;
    title_obs = nothing,
)
    # Build one plot for each variable of the current model.
    #
    # Arguments:
    #
    #     grid      - Makie GridLayout where the axes will be placed
    #     sim       - current simulation state
    #     title_obs - optional observable used as the title of the first axis
    #
    # Returns:
    #
    #     PlotPanel
    #
    # The plot panel stores:
    #
    #     axes        - one Axis per variable
    #     observables - one Observable per variable

    model = sim.model
    U = solution_matrix(sim)

    axes = Axis[]
    observables = Observable{Vector{Float64}}[]

    for j in 1:model.nvars
        title = if j == 1 && title_obs !== nothing
            title_obs
        else
            model.varnames[j]
        end

        ax = Axis(
            grid[j, 1],
            xlabel = j == model.nvars ? "x" : "",
            ylabel = model.varnames[j],
            title = title,
        )

        obs = Observable(copy(U[:, j]))

        lines!(ax, sim.x, obs; linewidth = 2)

        push!(axes, ax)
        push!(observables, obs)
    end

    panel = PlotPanel(axes, observables)

    rescale_solution_axes!(panel, sim)

    return panel
end


function refresh_plot_panel!(panel::PlotPanel, sim::SimulationState)
    # Refresh all solution plots from the current simulation state.

    model = sim.model
    U = solution_matrix(sim)

    length(panel.observables) == model.nvars ||
        error("Plot panel does not match the number of model variables.")

    for j in 1:model.nvars
        panel.observables[j][] = copy(U[:, j])
    end

    rescale_solution_axes!(panel, sim)

    return nothing
end


function rescale_solution_axes!(panel::PlotPanel, sim::SimulationState)
    # Automatically rescale each y-axis to the current range of its variable.

    model = sim.model
    U = solution_matrix(sim)

    length(panel.axes) == model.nvars ||
        error("Plot panel does not match the number of model variables.")

    for j in 1:model.nvars
        u = @view U[:, j]

        umin = minimum(u)
        umax = maximum(u)

        margin = max(0.2 * (umax - umin), 0.1)

        ylims!(
            panel.axes[j],
            umin - margin,
            umax + margin,
        )
    end

    return nothing
end


function relabel_plot_panel!(panel::PlotPanel, sim::SimulationState)
    # Update axis labels after changing the selected model.

    model = sim.model

    length(panel.axes) == model.nvars ||
        error("Plot panel does not match the number of model variables.")

    for j in 1:model.nvars
        panel.axes[j].ylabel = model.varnames[j]
        panel.axes[j].xlabel = j == model.nvars ? "x" : ""
    end

    return nothing
end