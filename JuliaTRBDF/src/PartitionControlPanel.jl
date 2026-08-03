# src/PartitionControlPanel.jl

# ============================================================
# Domain split / merge controls
# ============================================================


function update_split_marker!(
    app::AppState,
    segment::Int,
    left_count::Int,
)
    for marker in app.plot_panel.split_marker_observables
        marker[] = [NaN]
    end

    1 <= segment <= length(app.simulations) || return nothing
    sim = app.simulations[segment]
    2 <= left_count <= sim.N - 2 || return nothing

    displayed_length =
        segment_base_length(app, segment) *
        app.plot_panel.domain_length_scale
    split_position = displayed_length * left_count / sim.N
    app.plot_panel.split_marker_observables[segment][] = [split_position]

    return nothing
end


function rebuild_plot_panel_for_partition!(
    app::AppState,
    plot_grid::GridLayout;
    title_obs,
    domain_length_scale::Float64,
)
    previous_columns = length(app.plot_panel.segment_axes)
    clear_plot_panel!(app.plot_panel)

    app.plot_panel = build_plot_panel!(
        plot_grid,
        app;
        title_obs = title_obs,
    )

    set_plot_domain_scale!(
        app.plot_panel,
        app,
        domain_length_scale^2,
    )

    for column in (length(app.simulations) + 1):previous_columns
        colsize!(plot_grid, column, Fixed(0))
    end

    refresh_app_from_live_state!(app)

    return nothing
end


function split_domain_segment_app!(
    app::AppState,
    plot_grid::GridLayout,
    segment::Int,
    left_count::Int;
    title_obs,
)
    stop_worker!(app; wait = true)
    domain_length_scale = app.plot_panel.domain_length_scale

    lock(app.simlock)

    try
        app.generation[] += 1
        clear_snapshot_buffer!(app.snapshot_buffer)
        split_domain_segment!(app, segment, left_count)

        rebuild_plot_panel_for_partition!(
            app,
            plot_grid;
            title_obs = title_obs,
            domain_length_scale = domain_length_scale,
        )
    finally
        unlock(app.simlock)
    end

    return nothing
end


function merge_domain_segments_app!(
    app::AppState,
    plot_grid::GridLayout,
    left_segment::Int;
    title_obs,
)
    stop_worker!(app; wait = true)
    domain_length_scale = app.plot_panel.domain_length_scale

    lock(app.simlock)

    try
        app.generation[] += 1
        clear_snapshot_buffer!(app.snapshot_buffer)
        merge_domain_segments!(app, left_segment)

        rebuild_plot_panel_for_partition!(
            app,
            plot_grid;
            title_obs = title_obs,
            domain_length_scale = domain_length_scale,
        )
    finally
        unlock(app.simlock)
    end

    return nothing
end


function rebuild_partition_control_panel!(
    grid::GridLayout,
    app::AppState,
    item_ref,
    plot_grid::GridLayout;
    title_obs,
    selected_segment0::Int = 1,
)
    delete_control_items!(item_ref[])

    nsegments = length(app.simulations)
    selected_segment = Ref(clamp(selected_segment0, 1, nsegments))
    merge_rows = nsegments > 1 ? cld(nsegments - 1, 2) : 0
    last_control_row = nsegments > 1 ? 5 + merge_rows : 3

    title = Label(
        grid[1, 1:2],
        "Domain partition",
        tellwidth = false,
    )
    push!(item_ref[], title)

    layout_anchor = Label(
        grid[last_control_row, 2],
        "",
        tellwidth = false,
        tellheight = false,
        visible = false,
    )
    push!(item_ref[], layout_anchor)

    colsize!(grid, 1, Relative(0.5))
    colsize!(grid, 2, Relative(0.5))
    rowgap!(grid, 2)

    row_height = min(
        28.0,
        (212.0 - 2.0 * (last_control_row - 1)) / last_control_row,
    )

    for control_row in 1:last_control_row
        rowsize!(grid, control_row, Fixed(row_height))
    end

    row = 2
    segment_slider = nothing

    if nsegments > 1
        segment_label = Label(grid[row, 1], "Segment", tellwidth = false)
        segment_slider = Slider(
            grid[row, 2],
            range = 1:nsegments,
            startvalue = selected_segment[],
            tellwidth = false,
        )
        append!(item_ref[], Any[segment_label, segment_slider])
        row += 1
    end

    split_label_obs = Observable("Split point")
    split_label = Label(grid[row, 1], split_label_obs, tellwidth = false)
    selected_sim = app.simulations[selected_segment[]]
    can_split = selected_sim.N >= 4
    split_range = can_split ? (2:(selected_sim.N - 2)) : (1:1)
    initial_split = can_split ? clamp(div(selected_sim.N, 2), 2, selected_sim.N - 2) : 1
    split_slider = Slider(
        grid[row, 2],
        range = split_range,
        startvalue = initial_split,
        tellwidth = false,
    )
    append!(item_ref[], Any[split_label, split_slider])
    row += 1

    split_button = Button(
        grid[row, 1:2],
        label = "Split",
        tellwidth = false,
    )
    push!(item_ref[], split_button)
    row += 1

    function refresh_split_selection!()
        sim = app.simulations[selected_segment[]]
        can_split_now = sim.N >= 4

        if can_split_now
            split_slider.range[] = 2:(sim.N - 2)
            set_close_to!(split_slider, clamp(div(sim.N, 2), 2, sim.N - 2))
            left_count = Int(round(split_slider.value[]))
            split_label_obs[] = "Split point: $(left_count) / $(sim.N)"
            update_split_marker!(app, selected_segment[], left_count)
        else
            split_slider.range[] = 1:1
            split_label_obs[] = "Too few points"

            for marker in app.plot_panel.split_marker_observables
                marker[] = [NaN]
            end
        end

        return nothing
    end

    if segment_slider !== nothing
        on(segment_slider.value) do value
            selected_segment[] = Int(round(value))
            refresh_split_selection!()
            return nothing
        end
    end

    on(split_slider.value) do value
        sim = app.simulations[selected_segment[]]

        if sim.N >= 4
            left_count = Int(round(value))
            split_label_obs[] = "Split point: $(left_count) / $(sim.N)"
            update_split_marker!(app, selected_segment[], left_count)
        end

        return nothing
    end

    on(split_button.clicks) do _
        sim = app.simulations[selected_segment[]]
        sim.N >= 4 || return nothing
        left_count = Int(round(split_slider.value[]))
        split_at = selected_segment[]

        split_domain_segment_app!(
            app,
            plot_grid,
            split_at,
            left_count;
            title_obs = title_obs,
        )

        rebuild_partition_control_panel!(
            grid,
            app,
            item_ref,
            plot_grid;
            title_obs = title_obs,
            selected_segment0 = split_at,
        )

        return nothing
    end

    if nsegments > 1
        merge_title = Label(
            grid[row, 1:2],
            "Merge adjacent segments",
            tellwidth = false,
        )
        push!(item_ref[], merge_title)
        row += 1
        merge_start_row = row

        for boundary in 1:(nsegments - 1)
            left_boundary = boundary
            merge_row = merge_start_row + div(boundary - 1, 2)
            merge_column = 1 + mod(boundary - 1, 2)
            merge_button = Button(
                grid[merge_row, merge_column],
                label = "Merge $(boundary)–$(boundary + 1)",
                tellwidth = false,
            )
            push!(item_ref[], merge_button)

            on(merge_button.clicks) do _
                merge_domain_segments_app!(
                    app,
                    plot_grid,
                    left_boundary;
                    title_obs = title_obs,
                )

                rebuild_partition_control_panel!(
                    grid,
                    app,
                    item_ref,
                    plot_grid;
                    title_obs = title_obs,
                    selected_segment0 = min(left_boundary, length(app.simulations)),
                )

                return nothing
            end
        end
    end

    refresh_split_selection!()

    return nothing
end
