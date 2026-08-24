# src/TopMenu.jl

# ============================================================
# Model and boundary-condition controls
# ============================================================
#
# Selection area displayed at the top of the control window.
#
# Contains:
#
#     - model selector,
#     - boundary-condition selector.
#
# There are no separate text labels next to the dropdowns.
#
# ============================================================


function build_top_menu!(
    grid::GridLayout,
    app::AppState,
    plot_grid::GridLayout;
    registry,
    labels::Vector{String},
    first_label::String,
    N::Int,
    reltol::Float64,
    abstol::Float64,
    title_obs,
    model_name_obs::Observable{String},
    bc_name_obs::Observable{String},
)
    # Two dropdown menus placed side by side:
    #
    #     [model menu] [boundary condition menu]

    model_menu = Menu(
        grid[1, 1],
        options = labels,
        default = first_label,
        tellwidth = false,
    )

    bc_labels = boundary_condition_labels()

    boundary_menu = Menu(
        grid[1, 2],
        options = bc_labels,
        default = boundary_condition_label(app.sim.boundary_condition),
        tellwidth = false,
    )

    updating_boundary_menu = Ref(false)

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
            boundary_condition = app.sim.boundary_condition,
            title_obs = title_obs,
            model_name_obs = model_name_obs,
            bc_name_obs = bc_name_obs,
        )

        selected_boundary_label = boundary_condition_label(
            app.sim.boundary_condition,
        )

        if boundary_menu.selection[] != selected_boundary_label
            updating_boundary_menu[] = true

            try
                boundary_menu.selection[] = selected_boundary_label
            finally
                updating_boundary_menu[] = false
            end
        end
    end

    on(boundary_menu.selection) do selected_label
        updating_boundary_menu[] && return

        boundary_condition = boundary_condition_from_label(selected_label)

        switch_boundary_condition_app!(
            app,
            plot_grid,
            boundary_condition;
            N = N,
            dtmax = current_dtmax(app.sim),
            reltol = reltol,
            abstol = abstol,
            title_obs = title_obs,
            bc_name_obs = bc_name_obs,
        )
    end

    return (;
        model_menu,
        boundary_menu,
    )
end
