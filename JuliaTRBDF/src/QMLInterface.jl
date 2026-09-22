module ReactionDiffusionQML

get!(ENV, "QT_QUICK_CONTROLS_STYLE", "Basic")
get!(ENV, "QSG_RENDER_LOOP", "basic")

using GLMakie
import CairoMakie
using CxxWrap
using Base64
using LaTeXStrings
using Logging
using QML
using QMLMakie
using Printf

import ..ReactionDiffusionApp
const RD = ReactionDiffusionApp

const QML_FILE = normpath(joinpath(@__DIR__, "..", "qml", "Main.qml"))


mutable struct QMLBindings
    active_model_key::Observable{String}
    active_model_name::Observable{String}
    domain_length::Observable{Float64}
    random_mode::Observable{Bool}
    absolute_mode::Observable{Bool}
    selected_segment::Observable{Int}
    segment_count::Observable{Int}
    split_index::Observable{Int}
    split_maximum::Observable{Int}
    variables_json::Observable{String}
    equation_images_json::Observable{String}
    equation_preferred_width::Observable{Float64}
    perturbation_width::Observable{Float64}
    perturbation_height::Observable{Float64}
    checkpoint_available::Observable{Bool}
    message::Observable{String}
    series_running::Observable{Bool}
    series_completed_runs::Observable{Int}
    series_total_runs::Observable{Int}
    series_status::Observable{String}
    series_run_count::Observable{Int}
    series_maximum_time::Observable{Float64}
    series_maximum_steps::Observable{Int}
    series_check_interval::Observable{Float64}
    series_tolerance::Observable{Float64}
    series_required_checks::Observable{Int}
    series_dtmax::Observable{Float64}
    series_seed::Observable{String}
    series_head_variable::Observable{Int}
    series_live_preview::Observable{Bool}
    series_selected_segment::Observable{Int}
    series_selected_variable::Observable{Int}
    series_position::Observable{Float64}
    series_selected_panel_length::Observable{Float64}
    series_perturbations_json::Observable{String}
    series_results_json::Observable{String}
    series_mode::Observable{Bool}
end


mutable struct SeriesController
    settings::RD.SeriesSettings
    perturbations::Vector{RD.SeriesPerturbation}
    next_perturbation_id::Int
    selected_perturbation_id::Int
    selected_segment::Int
    selected_variable::Int
    selected_position::Float64
    editor_open::Bool
    preview_items::Vector{Any}
    running::Threads.Atomic{Bool}
    stop_requested::Threads.Atomic{Bool}
    task_ref::Base.RefValue{Union{Nothing, Task}}
    lock::ReentrantLock
    templates::Vector{RD.SeriesSegmentTemplate}
    base_snapshots::Vector{RD.SimulationSnapshot}
    latest_snapshots::Vector{Union{Nothing, RD.SimulationSnapshot}}
    snapshot_revision::Int
    published_snapshot_revision::Int
    generation::Int
    completed_runs::Int
    status::String
    location_counts::Vector{Vector{Int}}
    head_count_counts::Vector{Dict{Int, Int}}
    not_converged_counts::Vector{Int}
    results_revision::Int
    published_results_revision::Int
    finish_restores_base::Bool
end


mutable struct QMLController
    app::RD.AppState
    plot_grid::GridLayout
    title_obs
    model_name_obs::Observable{String}
    boundary_name_obs::Observable{String}
    registry::Dict{String, RD.ModelSpec}
    equation_images_by_model::Dict{String, Vector{String}}
    equation_widths_by_model::Dict{String, Float64}
    bindings::QMLBindings
    series::SeriesController
    diffusion_scale::Float64
    selected_segment::Int
    split_index::Int
    steps_per_frame::Int
    worker_sleep_time::Float64
    reltol::Float64
    abstol::Float64
    graphics_actions::Vector{Function}
    graphics_actions_lock::ReentrantLock
    graphics_busy::Observable{Bool}
    close_requested::Threads.Atomic{Bool}
end


const ACTIVE_QML_CONTROLLER = Ref{Union{Nothing, QMLController}}(nothing)
const QML_RENDER_CFUNCTION = Ref{Any}(nothing)


function report_startup_stage(label::AbstractString, started_ns::UInt64)
    elapsed = (time_ns() - started_ns) / 1.0e9
    @printf("[startup] %s %8.3f s\n", rpad(label, 36, '.'), elapsed)
    flush(stdout)
    return nothing
end


function set_if_changed!(observable::Observable, value)
    observable[] == value || (observable[] = value)
    return nothing
end


function json_string(value::AbstractString)
    io = IOBuffer()
    print(io, '"')

    for character in String(value)
        if character == '"'
            print(io, "\\\"")
        elseif character == '\\'
            print(io, "\\\\")
        elseif character == '\n'
            print(io, "\\n")
        elseif character == '\r'
            print(io, "\\r")
        elseif character == '\t'
            print(io, "\\t")
        elseif Int(character) < 0x20
            @printf(io, "\\u%04x", Int(character))
        else
            print(io, character)
        end
    end

    print(io, '"')
    return String(take!(io))
end


function json_string_array(values)
    return "[" * join(json_string.(String.(collect(values))), ",") * "]"
end



function json_number(value::Real)
    number = Float64(value)
    return isfinite(number) ? @sprintf("%.17g", number) : "null"
end


function empty_series_controller()
    return SeriesController(
        RD.SeriesSettings(),
        RD.SeriesPerturbation[],
        1,
        0,
        1,
        1,
        0.5,
        false,
        Any[],
        Threads.Atomic{Bool}(false),
        Threads.Atomic{Bool}(false),
        Ref{Union{Nothing, Task}}(nothing),
        ReentrantLock(),
        RD.SeriesSegmentTemplate[],
        RD.SimulationSnapshot[],
        Union{Nothing, RD.SimulationSnapshot}[],
        0,
        -1,
        0,
        0,
        "Configure perturbations",
        Vector{Int}[],
        Vector{Dict{Int, Int}}[],
        Int[],
        0,
        -1,
        false,
    )
end


function series_perturbations_json(series::SeriesController)
    entries = String[]

    for perturbation in series.perturbations
        push!(
            entries,
            "{" *
            "\"id\":" * string(perturbation.id) * "," *
            "\"panel\":" * string(perturbation.segment) * "," *
            "\"variable\":" * string(perturbation.variable) * "," *
            "\"position\":" * json_number(perturbation.position) * "," *
            "\"widthMin\":" * json_number(perturbation.width_min) * "," *
            "\"widthMax\":" * json_number(perturbation.width_max) * "," *
            "\"heightMin\":" * json_number(perturbation.height_min) * "," *
            "\"heightMax\":" * json_number(perturbation.height_max) *
            "}",
        )
    end

    return "[" * join(entries, ",") * "]"
end


function series_results_json(controller)
    series = controller.series
    app = controller.app
    entries = String[]

    lock(series.lock)

    try
        for segment in eachindex(series.location_counts)
            location_counts = series.location_counts[segment]
            count_map = series.head_count_counts[segment]
            maximum_heads = isempty(count_map) ? 0 : maximum(keys(count_map))
            head_labels = [json_string(string(count)) for count in 0:maximum_heads]
            head_values = [string(get(count_map, count, 0)) for count in 0:maximum_heads]
            push!(head_labels, json_string("Not converged"))
            push!(head_values, string(series.not_converged_counts[segment]))
            # Every partition panel displays local coordinates.  Keep the
            # histogram axis consistent with that convention, including the
            # right endpoint of a periodic panel (which is not itself a node).
            xmin = 0.0
            xmax = segment <= length(app.simulations) ?
                series_display_length(app, segment) :
                1.0

            push!(
                entries,
                "{" *
                "\"panel\":" * string(segment) * "," *
                "\"locationCounts\":[" * join(string.(location_counts), ",") * "]," *
                "\"xMin\":" * json_number(xmin) * "," *
                "\"xMax\":" * json_number(xmax) * "," *
                "\"headLabels\":[" * join(head_labels, ",") * "]," *
                "\"headCounts\":[" * join(head_values, ",") * "]" *
                "}",
            )
        end
    finally
        unlock(series.lock)
    end

    return "[" * join(entries, ",") * "]"
end


function series_perturbation_by_id(series::SeriesController, id::Int)
    index = findfirst(perturbation -> perturbation.id == id, series.perturbations)
    return index === nothing ? nothing : series.perturbations[index]
end


function series_display_length(app::RD.AppState, segment::Int)
    return RD.segment_base_length(app, segment) * app.plot_panel.domain_length_scale
end


function series_default_position(app::RD.AppState, segment::Int)
    return series_display_length(app, segment) / 2
end


const SERIES_WIDTH_DIGITS = 2
const SERIES_WIDTH_STEP = 10.0^-SERIES_WIDTH_DIGITS
const SERIES_HEIGHT_DIGITS = 1


function clamp_series_perturbation!(app::RD.AppState, perturbation::RD.SeriesPerturbation)
    1 <= perturbation.segment <= length(app.simulations) || return false
    displayed_length = series_display_length(app, perturbation.segment)
    perturbation.position = clamp(perturbation.position, 0.0, displayed_length)
    perturbation.width_max = clamp(perturbation.width_max, eps(Float64), displayed_length)
    perturbation.width_min = clamp(perturbation.width_min, eps(Float64), perturbation.width_max)
    # Keep the stored values at the precision shown in the series window, so
    # that focusing and leaving a field never changes a perturbation silently.
    perturbation.width_max = max(
        SERIES_WIDTH_STEP,
        round(perturbation.width_max; digits = SERIES_WIDTH_DIGITS),
    )
    perturbation.width_min = clamp(
        round(perturbation.width_min; digits = SERIES_WIDTH_DIGITS),
        SERIES_WIDTH_STEP,
        perturbation.width_max,
    )
    perturbation.height_min = round(perturbation.height_min; digits = SERIES_HEIGHT_DIGITS)
    perturbation.height_max = round(perturbation.height_max; digits = SERIES_HEIGHT_DIGITS)
    perturbation.variable = clamp(
        perturbation.variable,
        1,
        app.simulations[perturbation.segment].model.nvars,
    )
    return true
end

function render_model_equations_svg_uri(equations::AbstractVector{<:AbstractString})
    equation_lines = [
        filter(!isempty, strip.(split(String(equation), '\n')))
        for equation in equations
    ]
    filter!(!isempty, equation_lines)
    isempty(equation_lines) && return (uri = "", width = 0.0)

    figure = CairoMakie.Figure(
        figure_padding = (18, 18, 14, 14),
        backgroundcolor = :transparent,
    )
    row = 1

    for lines in equation_lines
        for (line_index, line) in enumerate(lines)
            parts = split(line, '='; limit = 2)
            has_equals = line_index == 1 && length(parts) == 2
            left_side = has_equals ? strip(parts[1]) : ""
            right_side = has_equals ? strip(parts[2]) : line

            if !isempty(left_side)
                CairoMakie.Label(
                    figure[row, 1],
                    LaTeXStrings.latexstring(left_side);
                    fontsize = 24,
                    color = :black,
                    halign = :right,
                    tellwidth = true,
                    tellheight = true,
                )
                CairoMakie.Label(
                    figure[row, 2],
                    LaTeXStrings.latexstring("=");
                    fontsize = 24,
                    color = :black,
                    halign = :center,
                    tellwidth = false,
                    tellheight = true,
                )
            end

            CairoMakie.Label(
                figure[row, 3],
                LaTeXStrings.latexstring(right_side);
                fontsize = 24,
                color = :black,
                halign = :left,
                tellwidth = true,
                tellheight = true,
            )
            row += 1
        end
    end

    colsize!(figure.layout, 1, Auto())
    colsize!(figure.layout, 2, Fixed(58))
    colsize!(figure.layout, 3, Auto())
    colgap!(figure.layout, 10)
    rowgap!(figure.layout, 18)
    CairoMakie.resize_to_layout!(figure)

    rendered_width = Float64(figure.scene.viewport[].widths[1])
    svg = sprint(show, MIME"image/svg+xml"(), figure)
    return (
        uri = "data:image/svg+xml;base64," * base64encode(svg),
        width = rendered_width,
    )
end


function render_equation_catalog(registry::Dict{String, RD.ModelSpec})
    images = Dict{String, Vector{String}}()
    widths = Dict{String, Float64}()
    CairoMakie.activate!(type = "svg")

    try
        for (key, model) in registry
            rendered = render_model_equations_svg_uri(model.latex_equations)
            images[key] = isempty(rendered.uri) ? String[] : String[rendered.uri]
            widths[key] = rendered.width
        end
    finally
        GLMakie.activate!()
    end

    return images, widths
end


function model_catalog_json(registry::Dict{String, RD.ModelSpec})
    families = String[]

    for family in RD.model_menu_catalog(registry)
        models = [
            "{" *
            "\"key\":" * json_string(entry.key) * "," *
            "\"label\":" * json_string(entry.label) *
            "}"
            for entry in family.models
        ]
        push!(
            families,
            "{" *
            "\"family\":" * json_string(family.family) * "," *
            "\"models\":[" * join(models, ",") * "]" *
            "}",
        )
    end

    return "[" * join(families, ",") * "]"
end


function active_model_menu_name(key::AbstractString)
    return "$(RD.model_family_name(key)) / $(RD.model_variant_name(key))"
end


function current_perturbation_modes(app::RD.AppState)
    state = RD.active_perturbation_state(app)
    state === nothing && return (true, false)

    return (state.random_mode[], state.absolute_mode[])
end


function update_model_bindings!(controller::QMLController)
    model = controller.app.sim.model
    bindings = controller.bindings
    active_key = bindings.active_model_key[]
    set_if_changed!(bindings.active_model_name, active_model_menu_name(active_key))
    set_if_changed!(bindings.variables_json, json_string_array(model.varnames))
    set_if_changed!(
        bindings.equation_images_json,
        json_string_array(get(controller.equation_images_by_model, active_key, String[])),
    )
    set_if_changed!(
        bindings.equation_preferred_width,
        get(controller.equation_widths_by_model, active_key, 0.0),
    )

    return nothing
end


function update_partition_bindings!(
    controller::QMLController;
    reset_index::Bool = false,
)
    app = controller.app
    count = length(app.simulations)
    controller.selected_segment = clamp(controller.selected_segment, 1, count)
    simulation = app.simulations[controller.selected_segment]
    split_maximum = max(2, simulation.N - 2)

    if reset_index || !(2 <= controller.split_index <= split_maximum)
        controller.split_index = clamp(round(Int, simulation.N / 2), 2, split_maximum)
    end

    set_if_changed!(controller.bindings.selected_segment, controller.selected_segment)
    set_if_changed!(controller.bindings.segment_count, count)
    set_if_changed!(controller.bindings.split_maximum, split_maximum)
    set_if_changed!(controller.bindings.split_index, controller.split_index)

    return nothing
end


function refresh_qml_state!(controller::QMLController)
    controller.graphics_busy[] && return nothing

    RD.refresh_ui_from_latest_snapshots!(controller.app)
    random_mode, absolute_mode = current_perturbation_modes(controller.app)
    perturbation_values = RD.perturbation_control_values(controller.app)
    set_if_changed!(controller.bindings.random_mode, random_mode)
    set_if_changed!(controller.bindings.absolute_mode, absolute_mode)
    set_if_changed!(
        controller.bindings.perturbation_width,
        perturbation_values.width,
    )
    set_if_changed!(
        controller.bindings.perturbation_height,
        perturbation_values.height,
    )
    set_if_changed!(
        controller.bindings.domain_length,
        controller.app.plot_panel.domain_length_scale,
    )
    set_if_changed!(
        controller.bindings.checkpoint_available,
        controller.app.saved_state[] !== nothing,
    )
    update_partition_bindings!(controller)

    return nothing
end


function enqueue_graphics_action!(
    controller::QMLController,
    action::Function,
)
    lock(controller.graphics_actions_lock)

    try
        controller.graphics_busy[] && return false
        controller.graphics_busy[] = true
        push!(controller.graphics_actions, action)
    finally
        unlock(controller.graphics_actions_lock)
    end

    return true
end


function take_graphics_action!(controller::QMLController)
    lock(controller.graphics_actions_lock)

    try
        isempty(controller.graphics_actions) && return nothing
        return popfirst!(controller.graphics_actions)
    finally
        unlock(controller.graphics_actions_lock)
    end
end


function process_graphics_actions!(controller::QMLController)
    while true
        action = take_graphics_action!(controller)
        action === nothing && break

        try
            action()
        catch error
            report_error!(controller, "Graphics update failed", error)
        finally
            controller.graphics_busy[] = false
        end

        refresh_qml_state!(controller)
    end

    return nothing
end


function qml_renderfunction(screen, scene_or_figure)
    controller = ACTIVE_QML_CONTROLLER[]

    if controller !== nothing
        process_graphics_actions!(controller)
    end

    QMLMakie.renderfunction(screen, scene_or_figure)
    return nothing
end


function install_qml_renderfunction!(controller::QMLController)
    ACTIVE_QML_CONTROLLER[] = controller

    if QML_RENDER_CFUNCTION[] === nothing
        QML_RENDER_CFUNCTION[] =
            @safe_cfunction(qml_renderfunction, Cvoid, (Any, Any))
    end

    QML.set_default_makie_renderfunction(QML_RENDER_CFUNCTION[])
    return nothing
end


function report_error!(controller::QMLController, context::AbstractString, error)
    message = "$context: $(sprint(showerror, error))"
    controller.bindings.message[] = message
    @error context exception = (error, catch_backtrace())

    return false
end


function guarded_action(action, controller::QMLController, context::AbstractString)
    try
        action()
        set_if_changed!(controller.bindings.message, "")
        refresh_qml_state!(controller)
        return true
    catch error
        return report_error!(controller, context, error)
    end
end


# ============================================================
# Series controls and runtime bridge
# ============================================================

function refresh_series_bindings!(controller::QMLController)
    series = controller.series
    settings = series.settings
    set_if_changed!(controller.bindings.series_running, series.running[])
    set_if_changed!(controller.bindings.series_completed_runs, series.completed_runs)
    set_if_changed!(controller.bindings.series_total_runs, settings.run_count)
    set_if_changed!(controller.bindings.series_status, series.status)
    set_if_changed!(controller.bindings.series_run_count, settings.run_count)
    set_if_changed!(controller.bindings.series_maximum_time, settings.maximum_time)
    set_if_changed!(
        controller.bindings.series_maximum_steps,
        settings.maximum_steps_per_panel,
    )
    set_if_changed!(controller.bindings.series_check_interval, settings.check_interval)
    set_if_changed!(controller.bindings.series_tolerance, settings.tolerance)
    set_if_changed!(
        controller.bindings.series_required_checks,
        settings.required_consecutive_checks,
    )
    set_if_changed!(controller.bindings.series_dtmax, settings.dtmax)
    set_if_changed!(controller.bindings.series_seed, string(settings.seed))
    set_if_changed!(controller.bindings.series_head_variable, settings.head_variable)
    set_if_changed!(controller.bindings.series_live_preview, settings.live_preview)
    set_if_changed!(controller.bindings.series_selected_segment, series.selected_segment)
    set_if_changed!(controller.bindings.series_selected_variable, series.selected_variable)
    set_if_changed!(controller.bindings.series_position, series.selected_position)
    set_if_changed!(
        controller.bindings.series_selected_panel_length,
        series_display_length(controller.app, series.selected_segment),
    )
    set_if_changed!(
        controller.bindings.series_perturbations_json,
        series_perturbations_json(series),
    )

    if series.published_results_revision != series.results_revision
        set_if_changed!(
            controller.bindings.series_results_json,
            series_results_json(controller),
        )
        series.published_results_revision = series.results_revision
    end

    return nothing
end


function refresh_series_runtime!(controller::QMLController)
    series = controller.series
    snapshots = RD.SimulationSnapshot[]
    running = series.running[]
    task = series.task_ref[]
    publish_snapshots = false
    commit_partial_state = false

    lock(series.lock)
    try
        if (series.settings.live_preview || !running) &&
           series.snapshot_revision != series.published_snapshot_revision &&
           length(series.latest_snapshots) == length(controller.app.simulations) &&
           all(snapshot -> snapshot !== nothing, series.latest_snapshots)
            snapshots = RD.SimulationSnapshot[
                snapshot::RD.SimulationSnapshot for snapshot in series.latest_snapshots
            ]
            series.published_snapshot_revision = series.snapshot_revision
            publish_snapshots = true
        end

        commit_partial_state = !running &&
                               task !== nothing &&
                               istaskdone(task) &&
                               series.stop_requested[] &&
                               !series.finish_restores_base &&
                               length(series.latest_snapshots) ==
                               length(controller.app.simulations) &&
                               all(snapshot -> snapshot !== nothing, series.latest_snapshots)

        if commit_partial_state && isempty(snapshots)
            snapshots = RD.SimulationSnapshot[
                snapshot::RD.SimulationSnapshot for snapshot in series.latest_snapshots
            ]
        end

        # This flag means that a terminal state has already been handled:
        # either the base state after a normal finish or the partial state
        # after a user stop.
        commit_partial_state && (series.finish_restores_base = true)
    finally
        unlock(series.lock)
    end

    commit_partial_state && commit_series_partial_state!(controller, snapshots)

    if publish_snapshots &&
       length(snapshots) == length(controller.app.simulations) &&
       !isempty(snapshots)
        try
            snapshot = RD.partition_snapshot_from_segments(snapshots, series.generation)
            RD.refresh_app_from_snapshot!(controller.app, snapshot)
        catch error
            report_error!(controller, "Series display refresh failed", error)
        end
    end

    refresh_series_bindings!(controller)

    if !running && task !== nothing && istaskdone(task) && controller.graphics_busy[]
        controller.graphics_busy[] = false

        if series.editor_open
            show_series_previews!(controller)
            update_series_position_marker!(controller)
        end

        refresh_qml_state!(controller)
    end

    return nothing
end


function commit_series_partial_state!(
    controller::QMLController,
    snapshots::Vector{RD.SimulationSnapshot},
)
    app = controller.app
    length(snapshots) == length(app.simulations) || return nothing

    lock(app.simlock)
    try
        for segment in eachindex(app.simulations)
            sim = app.simulations[segment]
            snapshot = snapshots[segment]
            runtime = app.segment_runtimes[segment]

            lock(runtime.lock)
            try
                RD.restart_after_manual_change!(
                    sim,
                    copy(snapshot.y);
                    reltol = controller.reltol,
                    abstol = controller.abstol,
                )
                sim.time_offset[] = snapshot.t
                sim.step_counter[] = snapshot.steps
            finally
                unlock(runtime.lock)
            end
        end

        RD.clear_snapshot_buffer!(app.snapshot_buffer)
    finally
        unlock(app.simlock)
    end

    return nothing
end


function clear_series_previews!(controller::QMLController)
    series = controller.series

    for item in series.preview_items
        try
            RD.delete_plot_panel_item!(item)
        catch error
            @debug "Series preview item was already deleted." exception = error
        end
    end

    empty!(series.preview_items)
    RD.rescale_solution_axes!(controller.app.plot_panel, controller.app.simulations)
    return nothing
end


function series_preview_anchor(
    controller::QMLController,
    perturbation::RD.SeriesPerturbation,
)
    panel = controller.app.plot_panel
    x = panel.segment_x_observables[perturbation.segment][]
    y = panel.segment_observables[perturbation.segment][perturbation.variable][]
    index = RD.nearest_grid_index(x, perturbation.position)
    return y[index]
end


function add_series_preview_box!(
    controller::QMLController,
    perturbation::RD.SeriesPerturbation,
    width::Float64,
    height::Float64;
    color,
    linewidth::Float64,
)
    panel = controller.app.plot_panel
    axis = panel.segment_axes[perturbation.segment][perturbation.variable]
    base = series_preview_anchor(controller, perturbation)
    left = perturbation.position - width / 2
    right = perturbation.position + width / 2
    top = base + height
    x = [left, right, right, left, left]
    y = [base, base, top, top, base]
    outline = lines!(axis, x, y; color = color, linestyle = :dash, linewidth = linewidth)
    push!(controller.series.preview_items, outline)
    return nothing
end


function show_series_previews!(controller::QMLController)
    series = controller.series
    clear_series_previews!(controller)

    for perturbation in series.perturbations
        clamp_series_perturbation!(controller.app, perturbation) || continue
        add_series_preview_box!(
            controller,
            perturbation,
            perturbation.width_max,
            perturbation.height_max;
            color = (:gray, 0.40),
            linewidth = 1.5,
        )
        add_series_preview_box!(
            controller,
            perturbation,
            perturbation.width_min,
            perturbation.height_min;
            color = (:gray, 0.80),
            linewidth = 2.0,
        )
    end

    rescale_series_preview_axes!(controller)

    return nothing
end


function rescale_series_preview_axes!(controller::QMLController)
    app = controller.app
    panel = app.plot_panel
    nvars = app.sim.model.nvars

    for variable in 1:nvars
        values = Float64[]
        axes = Axis[]

        for segment in eachindex(app.simulations)
            append!(values, RD.finite_values(panel.segment_observables[segment][variable][]))
            append!(
                values,
                RD.finite_values(panel.segment_preview_observables[segment][variable][]),
            )

            for perturbation in controller.series.perturbations
                perturbation.segment == segment || continue
                perturbation.variable == variable || continue
                base = series_preview_anchor(controller, perturbation)
                push!(values, base, base + perturbation.height_min, base + perturbation.height_max)
            end

            push!(axes, panel.segment_axes[segment][variable])
        end

        RD.set_axes_y_limits_from_values!(axes, values)
    end

    return nothing
end


function update_series_position_marker!(controller::QMLController)
    series = controller.series
    segment = series.selected_segment
    1 <= segment <= length(controller.app.simulations) || return nothing
    sim = controller.app.simulations[segment]
    displayed_length = series_display_length(controller.app, segment)
    displayed_length > 0.0 && sim.N >= 5 || return nothing
    index = clamp(
        round(Int, series.selected_position / displayed_length * sim.N),
        2,
        sim.N - 2,
    )

    # Reuse the established split-marker implementation verbatim.  It owns
    # the persistent Makie plot, token-based fade and lifecycle handling.
    RD.update_split_marker!(controller.app, segment, index)
    controller.app.plot_panel.split_marker_color_observables[segment][] = :gold

    return nothing
end


function clear_series_position_marker!(controller::QMLController)
    RD.update_split_marker!(controller.app, 0, 0)
    return nothing
end


function update_series_editor_selection!(controller::QMLController)
    series = controller.series
    perturbation = series_perturbation_by_id(series, series.selected_perturbation_id)

    if perturbation !== nothing
        series.selected_segment = perturbation.segment
        series.selected_variable = perturbation.variable
        series.selected_position = perturbation.position
    end

    series.selected_segment = clamp(series.selected_segment, 1, length(controller.app.simulations))
    series.selected_variable = clamp(
        series.selected_variable,
        1,
        controller.app.simulations[series.selected_segment].model.nvars,
    )
    series.selected_position = clamp(
        series.selected_position,
        0.0,
        series_display_length(controller.app, series.selected_segment),
    )
    return nothing
end


function open_series_editor!(controller::QMLController)
    series = controller.series
    series.running[] && return nothing
    RD.stop_worker!(controller.app; wait = true)
    series.editor_open = true
    update_series_editor_selection!(controller)
    series.status = isempty(series.perturbations) ?
        "Add at least one perturbation" :
        "Ready to run"
    controller.series.editor_open && show_series_previews!(controller)
    refresh_series_bindings!(controller)
    return nothing
end


function close_series_editor!(controller::QMLController)
    controller.series.editor_open = false
    clear_series_previews!(controller)
    clear_series_position_marker!(controller)
    return nothing
end


function set_series_mode!(controller::QMLController, enabled)
    return guarded_action(controller, "Series mode change failed") do
        app = controller.app
        series = controller.series
        enabled = Bool(enabled)
        enabled == controller.bindings.series_mode[] && return nothing

        if enabled
            app.mouse_perturbations_enabled[] = false
            RD.clear_perturbation_previews!(app.plot_panel)
            controller.bindings.series_mode[] = true
            open_series_editor!(controller)
        else
            # Leaving series mode is blocked until the running series is stopped.
            series.running[] && return nothing
            controller.bindings.series_mode[] = false
            close_series_editor!(controller)
            app.mouse_perturbations_enabled[] = true
        end

        return nothing
    end
end


function select_series_segment!(controller::QMLController, segment_value)
    series = controller.series
    series.running[] && return nothing
    series.selected_segment = clamp(Int(segment_value), 1, length(controller.app.simulations))
    series.selected_variable = clamp(
        series.selected_variable,
        1,
        controller.app.simulations[series.selected_segment].model.nvars,
    )
    series.selected_position = clamp(
        series.selected_position,
        0.0,
        series_display_length(controller.app, series.selected_segment),
    )
    update_series_position_marker!(controller)
    refresh_series_bindings!(controller)
    return nothing
end


function select_series_variable!(controller::QMLController, variable_value)
    series = controller.series
    series.running[] && return nothing
    nvars = controller.app.simulations[series.selected_segment].model.nvars
    series.selected_variable = clamp(Int(variable_value), 1, nvars)
    update_series_position_marker!(controller)
    refresh_series_bindings!(controller)
    return nothing
end


function set_series_position!(controller::QMLController, value)
    return guarded_action(controller, "Series-position change failed") do
        series = controller.series
        series.running[] && return nothing
        displayed_length = series_display_length(
            controller.app,
            series.selected_segment,
        )
        series.selected_position = clamp(
            parse_finite_qml_number(value, "Series position"),
            0.0,
            displayed_length,
        )

        # The slider only chooses the centre for the *next* perturbation.
        # Existing perturbations and their dashed preview boxes are never
        # moved or recreated while the user drags it.
        update_series_position_marker!(controller)
        refresh_series_bindings!(controller)
    end
end


function add_series_perturbation!(controller::QMLController)
    series = controller.series
    series.running[] && return nothing
    segment = series.selected_segment
    variable = series.selected_variable
    displayed_length = series_display_length(controller.app, segment)
    perturbation = RD.SeriesPerturbation(
        id = series.next_perturbation_id,
        segment = segment,
        variable = variable,
        position = series.selected_position,
        width_min = 0.05 * displayed_length,
        width_max = 0.10 * displayed_length,
        height_min = 1.0,
        height_max = 2.0,
    )
    clamp_series_perturbation!(controller.app, perturbation)
    push!(series.perturbations, perturbation)
    series.next_perturbation_id += 1
    series.selected_perturbation_id = perturbation.id
    series.status = "Ready to run"
    series.editor_open && show_series_previews!(controller)
    clear_series_position_marker!(controller)
    refresh_series_bindings!(controller)
    return nothing
end


function select_series_perturbation!(controller::QMLController, id_value)
    series = controller.series
    series.running[] && return nothing
    perturbation = series_perturbation_by_id(series, Int(id_value))
    perturbation === nothing && return nothing
    series.selected_perturbation_id = perturbation.id
    update_series_editor_selection!(controller)
    update_series_position_marker!(controller)
    refresh_series_bindings!(controller)
    return nothing
end


function delete_series_perturbation!(controller::QMLController, id_value)
    series = controller.series
    series.running[] && return nothing
    id = Int(id_value)
    filter!(perturbation -> perturbation.id != id, series.perturbations)
    series.selected_perturbation_id == id && (series.selected_perturbation_id = 0)
    series.status = isempty(series.perturbations) ? "Add at least one perturbation" : "Ready to run"
    series.editor_open && show_series_previews!(controller)
    refresh_series_bindings!(controller)
    return nothing
end


function update_series_perturbation!(controller::QMLController, id_value, field_value, value)
    series = controller.series
    series.running[] && return nothing
    perturbation = series_perturbation_by_id(series, Int(id_value))
    perturbation === nothing && return nothing
    field = String(field_value)

    if field == "panel"
        old_length = series_display_length(controller.app, perturbation.segment)
        perturbation.segment = clamp(Int(value), 1, length(controller.app.simulations))
        new_length = series_display_length(controller.app, perturbation.segment)
        scale = old_length > 0.0 ? new_length / old_length : 1.0
        perturbation.position *= scale
        perturbation.width_min *= scale
        perturbation.width_max *= scale
    elseif field == "variable"
        perturbation.variable = Int(value)
    elseif field == "position"
        perturbation.position = parse_finite_qml_number(value, "Position")
    elseif field == "widthMin"
        perturbation.width_min = parse_finite_qml_number(value, "Width minimum")
    elseif field == "widthMax"
        perturbation.width_max = parse_finite_qml_number(value, "Width maximum")
    elseif field == "heightMin"
        perturbation.height_min = parse_finite_qml_number(value, "Height minimum")
    elseif field == "heightMax"
        perturbation.height_max = parse_finite_qml_number(value, "Height maximum")
    else
        error("Unknown series-perturbation field: $field")
    end

    clamp_series_perturbation!(controller.app, perturbation)
    perturbation.width_min <= perturbation.width_max ||
        error("Width minimum cannot exceed width maximum.")
    perturbation.height_min <= perturbation.height_max ||
        error("Height minimum cannot exceed height maximum.")
    series.selected_perturbation_id = perturbation.id
    update_series_editor_selection!(controller)
    series.editor_open && show_series_previews!(controller)
    refresh_series_bindings!(controller)
    return nothing
end


function set_series_integer_setting!(controller::QMLController, field::Symbol, value)
    series = controller.series
    series.running[] && return nothing
    integer = Int(round(parse_finite_qml_number(value, String(field))))
    integer >= 1 || error("$(field) must be at least one.")

    if field == :run_count
        series.settings.run_count = integer
    elseif field == :maximum_steps
        series.settings.maximum_steps_per_panel = integer
    elseif field == :required_checks
        series.settings.required_consecutive_checks = integer
    elseif field == :head_variable
        nvars = controller.app.sim.model.nvars
        series.settings.head_variable = clamp(integer, 1, nvars)
    else
        error("Unknown series integer setting: $field")
    end

    refresh_series_bindings!(controller)
    return nothing
end


function set_series_float_setting!(controller::QMLController, field::Symbol, value)
    series = controller.series
    series.running[] && return nothing
    number = parse_finite_qml_number(value, String(field))
    number > 0.0 || error("$(field) must be positive.")

    if field == :maximum_time
        series.settings.maximum_time = number
    elseif field == :check_interval
        series.settings.check_interval = number
    elseif field == :dtmax
        series.settings.dtmax = number
    elseif field == :tolerance
        series.settings.tolerance = number
    else
        error("Unknown series floating-point setting: $field")
    end

    refresh_series_bindings!(controller)
    return nothing
end


function set_series_seed!(controller::QMLController, value)
    series = controller.series
    series.running[] && return nothing
    parsed = tryparse(UInt64, strip(String(value)))
    parsed === nothing && error("Random seed must be a non-negative integer.")
    series.settings.seed = parsed
    refresh_series_bindings!(controller)
    return nothing
end


function set_series_live_preview!(controller::QMLController, enabled)
    controller.series.running[] && return nothing
    controller.series.settings.live_preview = Bool(enabled)
    refresh_series_bindings!(controller)
    return nothing
end


function series_runtime_perturbations(
    controller::QMLController,
    templates::Vector{RD.SeriesSegmentTemplate},
)
    runtime_perturbations = RD.SeriesPerturbation[]

    for perturbation in controller.series.perturbations
        template = templates[perturbation.segment]
        display_length = series_display_length(controller.app, perturbation.segment)
        display_length > 0.0 || error("Selected panel has zero displayed length.")
        solver_length = template.boundary_condition == :periodic ?
            last(template.x) - first(template.x) + template.dx :
            max(last(template.x) - first(template.x), 0.0)
        scale = solver_length / display_length
        push!(
            runtime_perturbations,
            RD.SeriesPerturbation(
                id = perturbation.id,
                segment = perturbation.segment,
                variable = perturbation.variable,
                position = first(template.x) + perturbation.position * scale,
                width_min = perturbation.width_min * scale,
                width_max = perturbation.width_max * scale,
                height_min = perturbation.height_min,
                height_max = perturbation.height_max,
            ),
        )
    end

    return runtime_perturbations
end


function capture_series_base!(controller::QMLController)
    app = controller.app
    RD.stop_worker!(app; wait = true)
    templates = RD.SeriesSegmentTemplate[]
    snapshots = RD.SimulationSnapshot[]

    lock(app.simlock)
    try
        app.generation[] += 1
        RD.clear_snapshot_buffer!(app.snapshot_buffer)
        generation = app.generation[]

        for segment in eachindex(app.simulations)
            runtime = app.segment_runtimes[segment]
            lock(runtime.lock)
            try
                push!(templates, RD.make_series_template(app.simulations[segment]))
                push!(snapshots, RD.make_snapshot(app.simulations[segment], generation))
            finally
                unlock(runtime.lock)
            end
        end

        return templates, snapshots, generation
    finally
        unlock(app.simlock)
    end
end


function record_series_outcomes!(
    controller::QMLController,
    outcomes::Vector{RD.SeriesPanelOutcome},
    templates::Vector{RD.SeriesSegmentTemplate},
)
    series = controller.series

    lock(series.lock)
    try
        for segment in eachindex(outcomes)
            outcome = outcomes[segment]
            outcome.converged || begin
                series.not_converged_counts[segment] += 1
                continue
            end
            head_count = length(outcome.heads)
            counts = series.head_count_counts[segment]
            counts[head_count] = get(counts, head_count, 0) + 1

            for head in outcome.heads
                index = RD.nearest_grid_index(templates[segment].x, head.position)
                series.location_counts[segment][index] += 1
            end
        end

        series.completed_runs += 1
        series.status = "Run $(series.completed_runs)/$(series.settings.run_count) complete"
        series.results_revision += 1
    finally
        unlock(series.lock)
    end

    return nothing
end


function start_series!(controller::QMLController)
    series = controller.series
    series.running[] && return nothing
    isempty(series.perturbations) && error("Add at least one perturbation before starting a series.")
    RD.validate_series_settings(series.settings, controller.app.sim.model.nvars)

    for perturbation in series.perturbations
        clamp_series_perturbation!(controller.app, perturbation) ||
            error("A perturbation points to an unavailable panel.")
    end

    templates, base_snapshots, generation = capture_series_base!(controller)
    runtime_perturbations = series_runtime_perturbations(controller, templates)
    clear_series_previews!(controller)
    clear_series_position_marker!(controller)

    lock(series.lock)
    try
        series.templates = templates
        series.base_snapshots = base_snapshots
        series.latest_snapshots = Union{Nothing, RD.SimulationSnapshot}[base_snapshots...]
        series.snapshot_revision += 1
        series.published_snapshot_revision = -1
        series.generation = generation
        series.completed_runs = 0
        series.status = "Running 0/$(series.settings.run_count)"
        series.location_counts = [zeros(Int, length(template.x)) for template in templates]
        series.head_count_counts = [Dict{Int, Int}() for _ in templates]
        series.not_converged_counts = zeros(Int, length(templates))
        series.results_revision += 1
        series.finish_restores_base = false
        series.stop_requested[] = false
        series.running[] = true
    finally
        unlock(series.lock)
    end

    controller.graphics_busy[] = true
    refresh_series_bindings!(controller)
    settings = deepcopy(series.settings)
    perturbations = deepcopy(runtime_perturbations)

    series.task_ref[] = Threads.@spawn begin
        stopped = false

        try
            for run_index in 1:settings.run_count
                series.stop_requested[] && (stopped = true; break)
                outcomes = RD.run_series_realization!(
                    templates,
                    perturbations,
                    settings,
                    generation,
                    run_index;
                    cancelled = () -> series.stop_requested[],
                    on_snapshot = (segment, snapshot) -> begin
                        lock(series.lock)
                        try
                            series.latest_snapshots[segment] = snapshot
                            series.snapshot_revision += 1
                        finally
                            unlock(series.lock)
                        end
                    end,
                )
                series.stop_requested[] && (stopped = true; break)
                record_series_outcomes!(controller, outcomes, templates)
            end
        catch error
            lock(series.lock)
            try
                series.status = "Series failed: $(sprint(showerror, error))"
            finally
                unlock(series.lock)
            end
            @error "Series worker failed." exception = (error, catch_backtrace())
        finally
            lock(series.lock)
            try
                if stopped
                    series.status = "Stopped after $(series.completed_runs)/$(settings.run_count) runs"
                elseif series.completed_runs == settings.run_count
                    series.status = "Completed $(settings.run_count) runs"
                    series.latest_snapshots = Union{Nothing, RD.SimulationSnapshot}[series.base_snapshots...]
                    series.snapshot_revision += 1
                    series.finish_restores_base = true
                end
                series.running[] = false
                series.results_revision += 1
            finally
                unlock(series.lock)
            end
        end
    end

    return nothing
end


function stop_series!(controller::QMLController)
    series = controller.series
    series.running[] || return nothing
    series.stop_requested[] = true
    series.status = "Stopping after the current safe solver step..."
    refresh_series_bindings!(controller)
    return nothing
end


function clear_series_perturbations!(controller::QMLController)
    series = controller.series
    empty!(series.perturbations)
    series.selected_perturbation_id = 0
    series.selected_segment = 1
    series.selected_variable = 1
    series.selected_position = series_default_position(controller.app, 1)
    series.status = "Add at least one perturbation"
    clear_series_previews!(controller)
    clear_series_position_marker!(controller)
    refresh_series_bindings!(controller)
    return nothing
end


function clamp_all_series_perturbations!(controller::QMLController)
    series = controller.series
    filter!(perturbation -> clamp_series_perturbation!(controller.app, perturbation), series.perturbations)
    update_series_editor_selection!(controller)
    return nothing
end


function remap_series_after_split!(
    controller::QMLController,
    split_segment::Int,
    split_position::Float64,
)
    for perturbation in controller.series.perturbations
        if perturbation.segment > split_segment
            perturbation.segment += 1
        elseif perturbation.segment == split_segment && perturbation.position > split_position
            perturbation.segment += 1
            perturbation.position -= split_position
        end
    end

    clamp_all_series_perturbations!(controller)
    controller.series.editor_open && show_series_previews!(controller)
    refresh_series_bindings!(controller)
    return nothing
end


function remap_series_after_merge!(
    controller::QMLController,
    left_segment::Int,
    left_display_length::Float64,
)
    right_segment = left_segment + 1

    for perturbation in controller.series.perturbations
        if perturbation.segment == right_segment
            perturbation.segment = left_segment
            perturbation.position += left_display_length
        elseif perturbation.segment > right_segment
            perturbation.segment -= 1
        end
    end

    clamp_all_series_perturbations!(controller)
    controller.series.editor_open && show_series_previews!(controller)
    refresh_series_bindings!(controller)
    return nothing
end


function remap_series_after_swap!(controller::QMLController, left_segment::Int)
    right_segment = left_segment + 1

    for perturbation in controller.series.perturbations
        if perturbation.segment == left_segment
            perturbation.segment = right_segment
        elseif perturbation.segment == right_segment
            perturbation.segment = left_segment
        end
    end

    clamp_all_series_perturbations!(controller)
    controller.series.editor_open && show_series_previews!(controller)
    refresh_series_bindings!(controller)
    return nothing
end


function remap_series_after_delete!(controller::QMLController, deleted_segment::Int)
    filter!(perturbation -> perturbation.segment != deleted_segment, controller.series.perturbations)

    for perturbation in controller.series.perturbations
        perturbation.segment > deleted_segment && (perturbation.segment -= 1)
    end

    clamp_all_series_perturbations!(controller)
    controller.series.editor_open && show_series_previews!(controller)
    refresh_series_bindings!(controller)
    return nothing
end


function rescale_series_perturbations!(controller::QMLController, factor::Float64)
    for perturbation in controller.series.perturbations
        perturbation.position *= factor
        perturbation.width_min *= factor
        perturbation.width_max *= factor
    end

    controller.series.selected_position *= factor
    clamp_all_series_perturbations!(controller)
    controller.series.editor_open && show_series_previews!(controller)
    controller.series.editor_open && update_series_position_marker!(controller)
    refresh_series_bindings!(controller)
    return nothing
end

function toggle_running!(controller::QMLController)
    return guarded_action(controller, "Simulation control failed") do
        app = controller.app

        if app.worker_running[] || app.running[] || app.synchronization_running[]
            RD.stop_worker!(app; wait = true)
            RD.update_all_perturbation_previews!(app; stop_simulation = false)
        else
            RD.clear_perturbation_previews!(app.plot_panel)
            RD.start_worker!(
                app;
                steps_per_frame = controller.steps_per_frame,
                sleep_time = controller.worker_sleep_time,
            )
        end
    end
end


function reset_simulation!(controller::QMLController)
    return guarded_action(controller, "Reset failed") do
        RD.reset_initial_condition_app!(
            controller.app;
            plot_grid = controller.plot_grid,
            title_obs = controller.title_obs,
            steps_per_frame = controller.steps_per_frame,
            worker_sleep_time = controller.worker_sleep_time,
        )
        controller.selected_segment = 1
        update_partition_bindings!(controller; reset_index = true)
        clear_series_perturbations!(controller)
    end
end


function save_current_state!(controller::QMLController)
    return guarded_action(controller, "State save failed") do
        RD.save_simulation_state!(controller.app)
        controller.bindings.checkpoint_available[] = true
    end
end


function restore_saved_state!(controller::QMLController)
    return guarded_action(controller, "State restore failed") do
        RD.restore_saved_simulation_state_app!(
            controller.app;
            plot_grid = controller.plot_grid,
            title_obs = controller.title_obs,
            reltol = controller.reltol,
            abstol = controller.abstol,
        )
        controller.diffusion_scale =
            controller.app.plot_panel.domain_length_scale^2
        controller.boundary_name_obs[] = RD.boundary_condition_label(
            controller.app.initial_boundary_condition,
        )
        controller.selected_segment = 1
        update_partition_bindings!(controller; reset_index = true)
        clear_series_perturbations!(controller)
    end
end


function select_model!(controller::QMLController, key)
    key_string = String(key)

    return guarded_action(controller, "Model change failed") do
        model = RD.get_model(controller.registry, key_string)
        RD.switch_model_app!(
            controller.app,
            controller.plot_grid,
            model;
            N = controller.app.initial_N,
            dtmax = RD.current_dtmax(controller.app.sim),
            reltol = controller.reltol,
            abstol = controller.abstol,
            boundary_condition = controller.app.sim.boundary_condition,
            title_obs = controller.title_obs,
            model_name_obs = controller.model_name_obs,
            bc_name_obs = controller.boundary_name_obs,
        )
        RD.set_diffusion_scale_app!(
            controller.app,
            controller.diffusion_scale;
            steps_per_frame = controller.steps_per_frame,
            worker_sleep_time = controller.worker_sleep_time,
        )
        controller.bindings.active_model_key[] = key_string
        controller.selected_segment = 1
        update_model_bindings!(controller)
        update_partition_bindings!(controller; reset_index = true)
        clear_series_perturbations!(controller)
    end
end


function select_boundary_condition!(controller::QMLController, label)
    label_string = String(label)

    return guarded_action(controller, "Boundary-condition change failed") do
        boundary_condition = RD.boundary_condition_from_label(label_string)
        RD.switch_boundary_condition_app!(
            controller.app,
            controller.plot_grid,
            boundary_condition;
            N = controller.app.initial_N,
            dtmax = RD.current_dtmax(controller.app.sim),
            reltol = controller.reltol,
            abstol = controller.abstol,
            title_obs = controller.title_obs,
            bc_name_obs = controller.boundary_name_obs,
        )
        RD.set_diffusion_scale_app!(
            controller.app,
            controller.diffusion_scale;
            steps_per_frame = controller.steps_per_frame,
            worker_sleep_time = controller.worker_sleep_time,
        )
        controller.selected_segment = 1
        update_partition_bindings!(controller; reset_index = true)
        clear_series_perturbations!(controller)
    end
end


function set_dt_exponent!(controller::QMLController, exponent)
    return guarded_action(controller, "Maximum time-step change failed") do
        RD.set_dtmax_app!(controller.app, 10.0^Float64(exponent))
    end
end


function set_domain_exponent!(controller::QMLController, exponent)
    return guarded_action(controller, "Domain rescale failed") do
        previous_display_scale = controller.app.plot_panel.domain_length_scale
        controller.diffusion_scale = 10.0^Float64(exponent)
        RD.set_diffusion_scale_app!(
            controller.app,
            controller.diffusion_scale;
            steps_per_frame = controller.steps_per_frame,
            worker_sleep_time = controller.worker_sleep_time,
        )
        current_display_scale = controller.app.plot_panel.domain_length_scale
        previous_display_scale > 0.0 &&
            rescale_series_perturbations!(
                controller,
                current_display_scale / previous_display_scale,
            )
    end
end


function apply_constant_initial_condition!(
    controller::QMLController,
    zero_based_variable,
    value_text,
)
    return guarded_action(controller, "Steady-state change failed") do
        value = tryparse(Float64, String(value_text))
        value === nothing && error("Invalid numeric value: $(String(value_text))")
        variable = Int(zero_based_variable) + 1
        RD.set_single_constant_initial_condition_app!(
            controller.app;
            variable = variable,
            value = value,
            segment = controller.selected_segment,
            steps_per_frame = controller.steps_per_frame,
            worker_sleep_time = controller.worker_sleep_time,
        )
    end
end


function update_selected_segment!(
    controller::QMLController,
    segment::Integer;
    show_split_marker::Bool,
)
    count = length(controller.app.simulations)
    controller.selected_segment = clamp(Int(segment), 1, count)
    update_partition_bindings!(controller; reset_index = true)

    if show_split_marker
        RD.update_split_marker!(
            controller.app,
            controller.selected_segment,
            controller.split_index,
        )
    end

    return nothing
end


function select_segment!(controller::QMLController, one_based_segment)
    return guarded_action(controller, "Segment selection failed") do
        update_selected_segment!(
            controller,
            Int(one_based_segment);
            show_split_marker = false,
        )
    end
end


function select_split_segment!(controller::QMLController, one_based_segment)
    return guarded_action(controller, "Split-panel selection failed") do
        update_selected_segment!(
            controller,
            Int(one_based_segment);
            show_split_marker = true,
        )
    end
end


function change_selected_segment!(controller::QMLController, direction)
    return guarded_action(controller, "Segment selection failed") do
        update_selected_segment!(
            controller,
            controller.selected_segment + Int(direction);
            show_split_marker = true,
        )
    end
end


function set_split_index!(controller::QMLController, index)
    return guarded_action(controller, "Split-point change failed") do
        simulation = controller.app.simulations[controller.selected_segment]
        controller.split_index = clamp(Int(round(index)), 2, simulation.N - 2)
        update_partition_bindings!(controller)
        RD.update_split_marker!(
            controller.app,
            controller.selected_segment,
            controller.split_index,
        )
    end
end


function split_selected_segment!(controller::QMLController)
    return guarded_action(controller, "Domain split failed") do
        split_segment = controller.selected_segment
        old_display_length = series_display_length(controller.app, split_segment)
        old_point_count = controller.app.simulations[split_segment].N
        split_position = old_display_length * controller.split_index / old_point_count
        clear_series_position_marker!(controller)
        success = RD.split_domain_segment_app!(
            controller.app,
            controller.plot_grid,
            split_segment,
            controller.split_index;
            title_obs = controller.title_obs,
        )
        success || error("The selected split point is not valid.")
        update_partition_bindings!(controller; reset_index = true)
        remap_series_after_split!(controller, split_segment, split_position)
    end
end


function merge_boundary!(controller::QMLController, one_based_boundary)
    return guarded_action(controller, "Domain merge failed") do
        boundary = Int(one_based_boundary)
        left_display_length = series_display_length(controller.app, boundary)
        clear_series_position_marker!(controller)
        success = RD.merge_domain_segments_app!(
            controller.app,
            controller.plot_grid,
            boundary;
            title_obs = controller.title_obs,
        )
        success || error("The selected boundary is not valid.")
        controller.selected_segment = clamp(
            controller.selected_segment,
            1,
            length(controller.app.simulations),
        )
        update_partition_bindings!(controller; reset_index = true)
        remap_series_after_merge!(controller, boundary, left_display_length)
    end
end


function swap_boundary!(controller::QMLController, one_based_boundary)
    return guarded_action(controller, "Domain swap failed") do
        boundary = Int(one_based_boundary)
        selected_before_swap = controller.selected_segment
        clear_series_position_marker!(controller)
        success = RD.swap_adjacent_domain_segments_app!(
            controller.app,
            controller.plot_grid,
            boundary;
            title_obs = controller.title_obs,
            steps_per_frame = controller.steps_per_frame,
            worker_sleep_time = controller.worker_sleep_time,
        )
        success || error("The selected swap boundary is not valid.")
        controller.selected_segment = if selected_before_swap == boundary
            boundary + 1
        elseif selected_before_swap == boundary + 1
            boundary
        else
            selected_before_swap
        end
        update_partition_bindings!(controller; reset_index = true)
        remap_series_after_swap!(controller, boundary)
    end
end


function delete_selected_segment!(controller::QMLController)
    return guarded_action(controller, "Domain deletion failed") do
        deleted_segment = controller.selected_segment
        clear_series_position_marker!(controller)
        success = RD.delete_domain_segment_app!(
            controller.app,
            controller.plot_grid,
            deleted_segment;
            title_obs = controller.title_obs,
            steps_per_frame = controller.steps_per_frame,
            worker_sleep_time = controller.worker_sleep_time,
        )
        success || error("The last remaining panel cannot be deleted.")
        controller.selected_segment = clamp(
            deleted_segment,
            1,
            length(controller.app.simulations),
        )
        update_partition_bindings!(controller; reset_index = true)
        remap_series_after_delete!(controller, deleted_segment)
    end
end


function synchronize_domains!(controller::QMLController)
    return guarded_action(controller, "Synchronization failed") do
        RD.synchronize_domains_app!(controller.app)
    end
end


function toggle_random_mode!(controller::QMLController)
    return guarded_action(controller, "Perturbation-mode change failed") do
        RD.toggle_perturbation_random_mode!(controller.app)
    end
end


function toggle_absolute_mode!(controller::QMLController)
    return guarded_action(controller, "Perturbation-mode change failed") do
        RD.toggle_perturbation_absolute_mode!(controller.app)
    end
end


function parse_finite_qml_number(value, label::AbstractString)
    number = value isa Real ? Float64(value) : tryparse(Float64, String(value))
    number === nothing && error("$label must be a number.")
    isfinite(number) || error("$label must be finite.")

    return number
end


function set_perturbation_width!(controller::QMLController, value)
    return guarded_action(controller, "Perturbation width change failed") do
        width = parse_finite_qml_number(value, "Width")
        0.0 < width <= 1.0 || error("Width must be greater than 0 and at most 1.")
        RD.set_perturbation_width_value!(controller.app, width) ||
            error("Perturbation controls are not available.")
    end
end


function set_perturbation_height!(controller::QMLController, value)
    return guarded_action(controller, "Perturbation height change failed") do
        height = parse_finite_qml_number(value, "Height")
        RD.set_perturbation_height_value!(controller.app, height) ||
            error("Perturbation controls are not available.")
    end
end


function request_close!(controller::QMLController)
    controller.close_requested[] = true
    return nothing
end


function shutdown!(controller::QMLController)
    stop_series!(controller)
    series_task = controller.series.task_ref[]
    if series_task !== nothing && series_task !== current_task() && !istaskdone(series_task)
        wait(series_task)
    end
    clear_series_previews!(controller)
    clear_series_position_marker!(controller)
    RD.stop_worker!(controller.app; wait = true)

    return nothing
end


function cleanup_qml_runtime!()
    try
        with_logger(NullLogger()) do
            GLMakie.closeall()
        end
    catch error
        @debug "Makie resources will be released with the QML context." exception = error
    end

    try
        QML.quit(QML.get_qmlengine())
    catch error
        @debug "QML engine was already closed." exception = error
    end

    try
        QML.quit()
    catch error
        @debug "QML application was already closed." exception = error
    end

    try
        QML.cleanup()
        QML.process_events()
    catch error
        @debug "QML resources were already released." exception = error
    end

    return nothing
end


function run_qml_event_loop!(
    controller::QMLController;
    poll_interval::Float64 = 0.015,
)
    while !controller.close_requested[]
        QML.process_eventloop_updates()
        QML.process_events()
        refresh_series_runtime!(controller)
        controller.close_requested[] && break

        yield()
        sleep(poll_interval)
    end

    return nothing
end


function register_qml_functions!(controller::QMLController)
    QML.qmlfunction("refreshUI", () -> refresh_qml_state!(controller))
    QML.qmlfunction("toggleRunning", () -> toggle_running!(controller))
    QML.qmlfunction("saveCurrentState", () -> save_current_state!(controller))
    QML.qmlfunction(
        "restoreSavedState",
        () -> enqueue_graphics_action!(
            controller,
            () -> restore_saved_state!(controller),
        ),
    )
    QML.qmlfunction(
        "resetSimulation",
        () -> enqueue_graphics_action!(
            controller,
            () -> reset_simulation!(controller),
        ),
    )
    QML.qmlfunction(
        "selectModel",
        key -> begin
            key_string = String(key)
            enqueue_graphics_action!(
                controller,
                () -> select_model!(controller, key_string),
            )
        end,
    )
    QML.qmlfunction(
        "selectBoundaryCondition",
        label -> begin
            label_string = String(label)
            enqueue_graphics_action!(
                controller,
                () -> select_boundary_condition!(controller, label_string),
            )
        end,
    )
    QML.qmlfunction("setDtExponent", value -> set_dt_exponent!(controller, value))
    QML.qmlfunction(
        "setDomainExponent",
        value -> set_domain_exponent!(controller, value),
    )
    QML.qmlfunction(
        "applyConstantInitialCondition",
        (index, value) -> apply_constant_initial_condition!(controller, index, value),
    )
    QML.qmlfunction(
        "selectSegment",
        index -> select_segment!(controller, index),
    )
    QML.qmlfunction(
        "selectSplitSegment",
        index -> select_split_segment!(controller, index),
    )
    QML.qmlfunction(
        "changeSelectedSegment",
        direction -> change_selected_segment!(controller, direction),
    )
    QML.qmlfunction("setSplitIndex", index -> set_split_index!(controller, index))
    QML.qmlfunction(
        "splitSelectedSegment",
        () -> enqueue_graphics_action!(
            controller,
            () -> split_selected_segment!(controller),
        ),
    )
    QML.qmlfunction(
        "mergeBoundary",
        index -> begin
            boundary = Int(index)
            enqueue_graphics_action!(
                controller,
                () -> merge_boundary!(controller, boundary),
            )
        end,
    )
    QML.qmlfunction(
        "swapBoundary",
        index -> begin
            boundary = Int(index)
            enqueue_graphics_action!(
                controller,
                () -> swap_boundary!(controller, boundary),
            )
        end,
    )
    QML.qmlfunction(
        "deleteSelectedSegment",
        () -> enqueue_graphics_action!(
            controller,
            () -> delete_selected_segment!(controller),
        ),
    )
    QML.qmlfunction("synchronizeDomains", () -> synchronize_domains!(controller))
    QML.qmlfunction("toggleRandomMode", () -> toggle_random_mode!(controller))
    QML.qmlfunction("toggleAbsoluteMode", () -> toggle_absolute_mode!(controller))
    QML.qmlfunction(
        "setPerturbationWidth",
        value -> set_perturbation_width!(controller, value),
    )
    QML.qmlfunction(
        "setPerturbationHeight",
        value -> set_perturbation_height!(controller, value),
    )
    QML.qmlfunction("setSeriesMode", value -> set_series_mode!(controller, value))
    QML.qmlfunction("selectSeriesSegment", value -> select_series_segment!(controller, value))
    QML.qmlfunction("selectSeriesVariable", value -> select_series_variable!(controller, value))
    QML.qmlfunction("setSeriesPosition", value -> set_series_position!(controller, value))
    QML.qmlfunction("addSeriesPerturbation", () -> add_series_perturbation!(controller))
    QML.qmlfunction("selectSeriesPerturbation", value -> select_series_perturbation!(controller, value))
    QML.qmlfunction("deleteSeriesPerturbation", value -> delete_series_perturbation!(controller, value))
    QML.qmlfunction(
        "updateSeriesPerturbation",
        (id, field, value) -> update_series_perturbation!(controller, id, field, value),
    )
    QML.qmlfunction(
        "setSeriesRunCount",
        value -> set_series_integer_setting!(controller, :run_count, value),
    )
    QML.qmlfunction(
        "setSeriesMaximumTime",
        value -> set_series_float_setting!(controller, :maximum_time, value),
    )
    QML.qmlfunction(
        "setSeriesMaximumSteps",
        value -> set_series_integer_setting!(controller, :maximum_steps, value),
    )
    QML.qmlfunction(
        "setSeriesCheckInterval",
        value -> set_series_float_setting!(controller, :check_interval, value),
    )
    QML.qmlfunction(
        "setSeriesTolerance",
        value -> set_series_float_setting!(controller, :tolerance, value),
    )
    QML.qmlfunction(
        "setSeriesRequiredChecks",
        value -> set_series_integer_setting!(controller, :required_checks, value),
    )
    QML.qmlfunction(
        "setSeriesDtmax",
        value -> set_series_float_setting!(controller, :dtmax, value),
    )
    QML.qmlfunction("setSeriesSeed", value -> set_series_seed!(controller, value))
    QML.qmlfunction(
        "setSeriesHeadVariable",
        value -> set_series_integer_setting!(controller, :head_variable, value),
    )
    QML.qmlfunction(
        "setSeriesLivePreview",
        value -> set_series_live_preview!(controller, value),
    )
    QML.qmlfunction("startSeries", () -> start_series!(controller))
    QML.qmlfunction("stopSeries", () -> stop_series!(controller))
    QML.qmlfunction("requestClose", () -> request_close!(controller))

    return nothing
end


function qml_property_map(
    controller::QMLController,
    catalog_json::String;
    auto_close_ms::Int = 0,
)
    app = controller.app
    bindings = controller.bindings

    return QML.JuliaPropertyMap(
        "running" => app.running,
        "synchronizationStatus" => app.synchronization_status,
        "dtmax" => app.dtmax_obs,
        "modelName" => bindings.active_model_name,
        "activeModelKey" => bindings.active_model_key,
        "boundaryName" => controller.boundary_name_obs,
        "domainLength" => bindings.domain_length,
        "randomMode" => bindings.random_mode,
        "absoluteMode" => bindings.absolute_mode,
        "perturbationWidth" => bindings.perturbation_width,
        "perturbationHeight" => bindings.perturbation_height,
        "checkpointAvailable" => bindings.checkpoint_available,
        "selectedSegment" => bindings.selected_segment,
        "segmentCount" => bindings.segment_count,
        "splitIndex" => bindings.split_index,
        "splitMaximum" => bindings.split_maximum,
        "variablesJson" => bindings.variables_json,
        "equationImagesJson" => bindings.equation_images_json,
        "equationPreferredWidth" => bindings.equation_preferred_width,
        "modelCatalogJson" => Observable(catalog_json),
        "message" => bindings.message,
        "graphicsBusy" => controller.graphics_busy,
        "seriesRunning" => bindings.series_running,
        "seriesCompletedRuns" => bindings.series_completed_runs,
        "seriesTotalRuns" => bindings.series_total_runs,
        "seriesStatus" => bindings.series_status,
        "seriesRunCount" => bindings.series_run_count,
        "seriesMaximumTime" => bindings.series_maximum_time,
        "seriesMaximumSteps" => bindings.series_maximum_steps,
        "seriesCheckInterval" => bindings.series_check_interval,
        "seriesTolerance" => bindings.series_tolerance,
        "seriesRequiredChecks" => bindings.series_required_checks,
        "seriesDtmax" => bindings.series_dtmax,
        "seriesSeed" => bindings.series_seed,
        "seriesHeadVariable" => bindings.series_head_variable,
        "seriesLivePreview" => bindings.series_live_preview,
        "seriesSelectedSegment" => bindings.series_selected_segment,
        "seriesSelectedVariable" => bindings.series_selected_variable,
        "seriesPosition" => bindings.series_position,
        "seriesSelectedPanelLength" => bindings.series_selected_panel_length,
        "seriesPerturbationsJson" => bindings.series_perturbations_json,
        "seriesResultsJson" => bindings.series_results_json,
        "seriesMode" => bindings.series_mode,
        "autoCloseMs" => Observable(auto_close_ms),
    )
end


function create_qml_controller(;
    N::Int,
    boundary_condition0::Symbol,
    dtmax0::Float64,
    reltol::Float64,
    abstol::Float64,
    steps_per_frame::Int,
    worker_sleep_time::Float64,
)
    RD.validate_boundary_condition(boundary_condition0)
    registry = RD.MODEL_REGISTRY
    stage_started_ns = time_ns()
    equation_images_by_model, equation_widths_by_model =
        render_equation_catalog(registry)
    report_startup_stage("Render equation catalog", stage_started_ns)

    labels = RD.model_labels(registry)
    isempty(labels) && error("No models found in MODEL_REGISTRY.")
    first_label = first(labels)
    first_key = RD.model_registry_key_for_label(registry, first_label)
    first_model = RD.get_model(registry, first_key)

    stage_started_ns = time_ns()
    simulation = RD.create_simulation_state(
        first_model;
        N = N,
        dtmax = dtmax0,
        reltol = reltol,
        abstol = abstol,
        boundary_condition = boundary_condition0,
    )
    report_startup_stage("Create initial simulation", stage_started_ns)

    stage_started_ns = time_ns()
    running = Observable(false)
    dtmax = Observable(dtmax0)
    dt = Observable(RD.current_internal_dt(simulation))
    time = Observable(RD.current_display_time(simulation))
    steps = Observable(simulation.step_counter[])
    model_name = Observable(first_model.display_name)
    boundary_name = Observable(RD.boundary_condition_label(boundary_condition0))
    title = lift(model_name, boundary_name, time, steps) do _, _, t, step_count
        "t = $(@sprintf("%.1e", t)) | steps = $step_count"
    end
    figure = Figure(size = (1300, 820))
    plot_grid = GridLayout(tellwidth = false, tellheight = false)
    figure[1, 1] = plot_grid
    rowsize!(figure.layout, 1, Auto(false, 1.0))
    colsize!(figure.layout, 1, Relative(1.0))
    app = RD.AppState(
        simulation,
        RD.SimulationState[simulation],
        N,
        boundary_condition0,
        RD.empty_plot_panel(),
        running,
        dtmax,
        dt,
        time,
        steps,
        Threads.Atomic{Bool}(false),
        Ref{Union{Nothing, Task}}(nothing),
        Ref{Union{Nothing, Task}}(nothing),
        RD.empty_snapshot_buffer(),
        Threads.Atomic{Int}(0),
        ReentrantLock(),
        RD.SegmentRuntime[RD.empty_segment_runtime()],
        Threads.Atomic{Bool}(false),
        Ref{Union{Nothing, Task}}(nothing),
        Observable("Synchronized"),
        Ref{Union{Nothing, RD.SavedSimulationState}}(nothing),
        false,
        Threads.Atomic{Bool}(true),
    )
    app.plot_panel = RD.build_plot_panel!(plot_grid, app; title_obs = title)
    report_startup_stage("Create application and plots", stage_started_ns)

    stage_started_ns = time_ns()
    series_defaults = RD.SeriesSettings()
    bindings = QMLBindings(
        Observable(first_key),
        Observable(active_model_menu_name(first_key)),
        Observable(1.0),
        Observable(true),
        Observable(false),
        Observable(1),
        Observable(1),
        Observable(clamp(round(Int, N / 2), 2, N - 2)),
        Observable(max(2, N - 2)),
        Observable(json_string_array(first_model.varnames)),
        Observable(json_string_array(equation_images_by_model[first_key])),
        Observable(get(equation_widths_by_model, first_key, 0.0)),
        Observable(0.05),
        Observable(0.0),
        Observable(false),
        Observable(""),
        Observable(false),
        Observable(0),
        Observable(100),
        Observable("Configure perturbations"),
        Observable(series_defaults.run_count),
        Observable(series_defaults.maximum_time),
        Observable(series_defaults.maximum_steps_per_panel),
        Observable(series_defaults.check_interval),
        Observable(series_defaults.tolerance),
        Observable(series_defaults.required_consecutive_checks),
        Observable(series_defaults.dtmax),
        Observable(string(series_defaults.seed)),
        Observable(1),
        Observable(false),
        Observable(1),
        Observable(1),
        Observable(0.5),
        Observable(1.0),
        Observable("[]"),
        Observable("[]"),
        Observable(false),
    )
    controller = QMLController(
        app,
        plot_grid,
        title,
        model_name,
        boundary_name,
        registry,
        equation_images_by_model,
        equation_widths_by_model,
        bindings,
        empty_series_controller(),
        1.0,
        1,
        bindings.split_index[],
        steps_per_frame,
        worker_sleep_time,
        reltol,
        abstol,
        Function[],
        ReentrantLock(),
        Observable(false),
        Threads.Atomic{Bool}(false),
    )
    controller.series.selected_position = series_default_position(controller.app, 1)
    refresh_series_bindings!(controller)
    report_startup_stage("Create QML controller", stage_started_ns)

    return controller, figure
end


function run_qml_app(;
    N::Int = 500,
    boundary_condition0::Symbol = :neumann,
    dtmax0::Float64 = 1e-2,
    reltol::Float64 = 1e-5,
    abstol::Float64 = 1e-7,
    steps_per_frame::Int = 5,
    worker_sleep_time::Float64 = 0.001,
    auto_close_ms::Int = 0,
    startup_started_ns::UInt64 = time_ns(),
)
    isfile(QML_FILE) || error("QML interface file does not exist: $QML_FILE")

    if Threads.nthreads() == 1
        @warn "Julia is running with one thread; use --threads=auto for independent solvers."
    end

    controller, figure = create_qml_controller(
        N = N,
        boundary_condition0 = boundary_condition0,
        dtmax0 = dtmax0,
        reltol = reltol,
        abstol = abstol,
        steps_per_frame = steps_per_frame,
        worker_sleep_time = worker_sleep_time,
    )
    stage_started_ns = time_ns()
    install_qml_renderfunction!(controller)
    report_startup_stage("Install QML render function", stage_started_ns)

    stage_started_ns = time_ns()
    register_qml_functions!(controller)
    report_startup_stage("Register QML callbacks", stage_started_ns)

    stage_started_ns = time_ns()
    properties = qml_property_map(
        controller,
        model_catalog_json(controller.registry),
        auto_close_ms = auto_close_ms,
    )
    report_startup_stage("Prepare QML properties", stage_started_ns)

    stage_started_ns = time_ns()
    QML.loadqml(QML_FILE; plot = figure, ui = properties)
    report_startup_stage("QML.loadqml", stage_started_ns)
    report_startup_stage("Total before QML event loop", startup_started_ns)
    println("[startup] Entering QML event loop")
    flush(stdout)

    try
        run_qml_event_loop!(controller)
    finally
        shutdown!(controller)
        cleanup_qml_runtime!()
        ACTIVE_QML_CONTROLLER[] = nothing
    end

    return nothing
end


export run_qml_app

end
