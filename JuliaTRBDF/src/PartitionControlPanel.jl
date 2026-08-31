# src/PartitionControlPanel.jl

# ============================================================
# Domain split / merge controls
# ============================================================

const PARTITION_WORKER_SETTLE_TIME = 0.05


function stop_and_settle_partition_workers!(app::AppState)
    stop_worker!(app; wait = true)

    # Give the UI/render tasks one short scheduling window after every solver
    # task has finished before replacing simulations and plot observables.
    yield()
    sleep(PARTITION_WORKER_SETTLE_TIME)

    return nothing
end


function update_split_marker!(
    app::AppState,
    segment::Int,
    left_count::Int,
)
    for index in eachindex(app.plot_panel.split_marker_observables)
        app.plot_panel.split_marker_fade_tokens[index][] += 1
        app.plot_panel.split_marker_observables[index][] = [NaN]
        app.plot_panel.split_marker_alpha_observables[index][] = 0.0
    end

    1 <= segment <= length(app.simulations) || return nothing
    sim = app.simulations[segment]
    2 <= left_count <= sim.N - 2 || return nothing

    displayed_length =
        segment_base_length(app, segment) *
        app.plot_panel.domain_length_scale
    split_position = displayed_length * left_count / sim.N
    marker = app.plot_panel.split_marker_observables[segment]
    alpha = app.plot_panel.split_marker_alpha_observables[segment]
    token_ref = app.plot_panel.split_marker_fade_tokens[segment]
    marker[] = [split_position]
    alpha[] = 1.0
    token = token_ref[]

    @async begin
        sleep(5.0)
        fade_duration = 0.5
        fade_started_at = time_ns()

        while true
            token_ref[] == token || return nothing
            elapsed = (time_ns() - fade_started_at) / 1e9
            elapsed >= fade_duration && break
            alpha[] = max(0.0, 1.0 - elapsed / fade_duration)
            sleep(1 / 60)
        end

        token_ref[] == token || return nothing
        marker[] = [NaN]
        alpha[] = 0.0
        return nothing
    end

    return nothing
end


function rebuild_plot_panel_for_partition!(
    app::AppState,
    plot_grid::GridLayout;
    title_obs,
    domain_length_scale::Float64,
)
    clear_plot_panel!(app.plot_panel)
    reset_plot_grid_layout!(plot_grid)

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
    lock(app.simlock)

    try
        1 <= segment <= length(app.simulations) || return false
        2 <= left_count <= app.simulations[segment].N - 2 || return false

        # Increment before the first yielding operation. Every callback from
        # the panel being replaced can then recognize that it is stale.
        app.generation[] += 1
        stop_and_settle_partition_workers!(app)
        domain_length_scale = app.plot_panel.domain_length_scale
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

    return true
end


function merge_domain_segments_app!(
    app::AppState,
    plot_grid::GridLayout,
    left_segment::Int;
    title_obs,
)
    lock(app.simlock)

    try
        1 <= left_segment < length(app.simulations) || return false

        # Invalidate the current control panel before waiting for workers.
        # This prevents a queued second click on the old Merge button from
        # starting another merge with obsolete segment indices.
        app.generation[] += 1
        stop_and_settle_partition_workers!(app)
        domain_length_scale = app.plot_panel.domain_length_scale

        synchronize_segment_indices_blocking!(
            app,
            [left_segment, left_segment + 1];
            allow_cancel = false,
        )

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

    return true
end


function finalize_partition_topology_change!(
    app::AppState,
    plot_grid::GridLayout;
    title_obs,
    domain_length_scale::Float64,
)
    rebuild_plot_panel_for_partition!(
        app,
        plot_grid;
        title_obs = title_obs,
        domain_length_scale = domain_length_scale,
    )

    snapshot = make_partition_snapshot(app.simulations, app.generation[])
    store_runtime_snapshots!(app, snapshot.segments)
    refresh_app_from_snapshot!(app, snapshot)

    return nothing
end


function swap_adjacent_domain_segments_app!(
    app::AppState,
    plot_grid::GridLayout,
    left_segment::Int;
    title_obs,
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
)
    1 <= left_segment < length(app.simulations) || return false
    was_running = app.worker_running[] ||
        any(runtime.running[] for runtime in app.segment_runtimes)

    lock(app.simlock)

    try
        app.generation[] += 1
        stop_and_settle_partition_workers!(app)
        domain_length_scale = app.plot_panel.domain_length_scale
        clear_snapshot_buffer!(app.snapshot_buffer)
        swap_adjacent_domain_segments!(app, left_segment)
        finalize_partition_topology_change!(
            app,
            plot_grid;
            title_obs = title_obs,
            domain_length_scale = domain_length_scale,
        )
    finally
        unlock(app.simlock)
    end

    if was_running
        start_worker!(
            app;
            steps_per_frame = steps_per_frame,
            sleep_time = worker_sleep_time,
        )
    end

    return true
end


function delete_domain_segment_app!(
    app::AppState,
    plot_grid::GridLayout,
    segment::Int;
    title_obs,
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
)
    length(app.simulations) > 1 || return false
    1 <= segment <= length(app.simulations) || return false
    was_running = app.worker_running[] ||
        any(runtime.running[] for runtime in app.segment_runtimes)

    lock(app.simlock)

    try
        app.generation[] += 1
        stop_and_settle_partition_workers!(app)
        domain_length_scale = app.plot_panel.domain_length_scale
        clear_snapshot_buffer!(app.snapshot_buffer)
        delete_domain_segment!(app, segment)
        finalize_partition_topology_change!(
            app,
            plot_grid;
            title_obs = title_obs,
            domain_length_scale = domain_length_scale,
        )
    finally
        unlock(app.simlock)
    end

    if was_running
        start_worker!(
            app;
            steps_per_frame = steps_per_frame,
            sleep_time = worker_sleep_time,
        )
    end

    return true
end


function rebuild_partition_control_panel!(
    grid::GridLayout,
    app::AppState,
    item_ref,
    plot_grid::GridLayout;
    title_obs,
    selected_segment0::Int = 1,
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
)
    delete_control_items!(item_ref[])

    nsegments = length(app.simulations)
    panel_generation = app.generation[]
    selected_segment = Ref(clamp(selected_segment0, 1, nsegments))
    updating_split_slider = Ref(false)
    selected_segment_label_obs = Observable(
        "Plot $(selected_segment[]) / $(nsegments)",
    )
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
    left_segment_button = nothing
    right_segment_button = nothing

    if nsegments > 1
        segment_selector_grid = GridLayout(
            tellwidth = true,
            tellheight = false,
        )
        grid[row, 1:2] = segment_selector_grid
        left_segment_button = Button(
            segment_selector_grid[1, 1],
            label = "Left",
            tellwidth = false,
        )
        selected_segment_label = Label(
            segment_selector_grid[1, 2],
            selected_segment_label_obs,
            tellwidth = false,
            halign = :center,
        )
        right_segment_button = Button(
            segment_selector_grid[1, 3],
            label = "Right",
            tellwidth = false,
        )
        colsize!(segment_selector_grid, 1, Fixed(72))
        colsize!(segment_selector_grid, 2, Fixed(105))
        colsize!(segment_selector_grid, 3, Fixed(72))
        colgap!(segment_selector_grid, 4)
        try
            segment_selector_grid.halign = :center
        catch
        end
        append!(
            item_ref[],
            Any[
                segment_selector_grid,
                left_segment_button,
                selected_segment_label,
                right_segment_button,
            ],
        )
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

    function refresh_split_selection!(; show_marker::Bool = false)
        sim = app.simulations[selected_segment[]]
        can_split_now = sim.N >= 4
        selected_segment_label_obs[] =
            "Plot $(selected_segment[]) / $(length(app.simulations))"

        for index in eachindex(app.plot_panel.split_marker_observables)
            app.plot_panel.split_marker_fade_tokens[index][] += 1
            app.plot_panel.split_marker_observables[index][] = [NaN]
            app.plot_panel.split_marker_alpha_observables[index][] = 0.0
        end

        if can_split_now
            updating_split_slider[] = true

            try
                split_slider.range[] = 2:(sim.N - 2)
                set_close_to!(split_slider, clamp(div(sim.N, 2), 2, sim.N - 2))
            finally
                updating_split_slider[] = false
            end

            left_count = Int(round(split_slider.value[]))
            split_label_obs[] = "Split point: $(left_count) / $(sim.N)"

            if show_marker
                update_split_marker!(app, selected_segment[], left_count)
            end
        else
            split_slider.range[] = 1:1
            split_label_obs[] = "Too few points"

            for marker in app.plot_panel.split_marker_observables
                marker[] = [NaN]
            end

            for alpha in app.plot_panel.split_marker_alpha_observables
                alpha[] = 0.0
            end
        end

        return nothing
    end

    if left_segment_button !== nothing
        on(left_segment_button.clicks) do _
            panel_generation == app.generation[] || return nothing
            selected_segment[] = max(1, selected_segment[] - 1)
            refresh_split_selection!(show_marker = true)
            return nothing
        end

        on(right_segment_button.clicks) do _
            panel_generation == app.generation[] || return nothing
            selected_segment[] = min(nsegments, selected_segment[] + 1)
            refresh_split_selection!(show_marker = true)
            return nothing
        end
    end

    on(split_slider.value) do value
        panel_generation == app.generation[] || return nothing
        updating_split_slider[] && return nothing
        sim = app.simulations[selected_segment[]]

        if sim.N >= 4
            left_count = Int(round(value))
            split_label_obs[] = "Split point: $(left_count) / $(sim.N)"
            update_split_marker!(app, selected_segment[], left_count)
        end

        return nothing
    end

    on(split_button.clicks) do _
        panel_generation == app.generation[] || return nothing
        sim = app.simulations[selected_segment[]]
        sim.N >= 4 || return nothing
        left_count = Int(round(split_slider.value[]))
        split_at = selected_segment[]

        did_split = split_domain_segment_app!(
            app,
            plot_grid,
            split_at,
            left_count;
            title_obs = title_obs,
        )
        did_split || return nothing

        rebuild_partition_control_panel!(
            grid,
            app,
            item_ref,
            plot_grid;
            title_obs = title_obs,
            selected_segment0 = split_at,
            steps_per_frame = steps_per_frame,
            worker_sleep_time = worker_sleep_time,
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
                panel_generation == app.generation[] || return nothing

                did_merge = merge_domain_segments_app!(
                    app,
                    plot_grid,
                    left_boundary;
                    title_obs = title_obs,
                )
                did_merge || return nothing

                rebuild_partition_control_panel!(
                    grid,
                    app,
                    item_ref,
                    plot_grid;
                    title_obs = title_obs,
                    selected_segment0 = min(left_boundary, length(app.simulations)),
                    steps_per_frame = steps_per_frame,
                    worker_sleep_time = worker_sleep_time,
                )

                return nothing
            end
        end
    end

    refresh_split_selection!()

    return nothing
end
