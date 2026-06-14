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
# In the threaded version, the UI should preferably update plots
# from SimulationSnapshot objects, not directly from the live
# integrator.
#
# ============================================================


function empty_plot_panel()
    # Create an empty plot panel.

    return PlotPanel(
        Axis[],
        Observable{Vector{Float64}}[],
    )
end


function solution_matrix_from_snapshot(snapshot::SimulationSnapshot)
    # Return the snapshot state reshaped as an N × nvars matrix.

    return reshape(snapshot.y, snapshot.N, snapshot.nvars)
end


function build_plot_panel!(
    grid::GridLayout,
    sim::SimulationState;
    title_obs = nothing,
)
    # Build one plot for each variable of the current model.
    #
    # This function is used when a model is initialized or switched.
    # At this moment the UI owns the operation, so reading directly
    # from sim is acceptable.

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
    #
    # This should only be used when the UI safely owns the simulation state,
    # for example immediately after initialization or after a locked manual
    # perturbation.

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


function refresh_plot_panel_from_snapshot!(
    panel::PlotPanel,
    snapshot::SimulationSnapshot,
)
    # Refresh all solution plots from a thread-safe simulation snapshot.
    #
    # This is the preferred update path while the worker thread is running.

    U = solution_matrix_from_snapshot(snapshot)

    length(panel.observables) == snapshot.nvars ||
        error("Plot panel does not match the number of snapshot variables.")

    for j in 1:snapshot.nvars
        panel.observables[j][] = copy(U[:, j])
    end

    rescale_solution_axes_from_snapshot!(panel, snapshot)

    return nothing
end


function rescale_solution_axes!(panel::PlotPanel, sim::SimulationState)
    # Automatically rescale each y-axis to the current range of its variable.
    #
    # This reads from the live simulation state, so it should not be called
    # concurrently with the worker thread.

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


function rescale_solution_axes_from_snapshot!(
    panel::PlotPanel,
    snapshot::SimulationSnapshot,
)
    # Automatically rescale each y-axis using a thread-safe snapshot.

    U = solution_matrix_from_snapshot(snapshot)

    length(panel.axes) == snapshot.nvars ||
        error("Plot panel does not match the number of snapshot variables.")

    for j in 1:snapshot.nvars
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