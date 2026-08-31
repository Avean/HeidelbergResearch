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
    end
end


function set_dt_exponent!(controller::QMLController, exponent)
    return guarded_action(controller, "Maximum time-step change failed") do
        RD.set_dtmax_app!(controller.app, 10.0^Float64(exponent))
    end
end


function set_domain_exponent!(controller::QMLController, exponent)
    return guarded_action(controller, "Domain rescale failed") do
        controller.diffusion_scale = 10.0^Float64(exponent)
        RD.set_diffusion_scale_app!(
            controller.app,
            controller.diffusion_scale;
            steps_per_frame = controller.steps_per_frame,
            worker_sleep_time = controller.worker_sleep_time,
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
        success = RD.split_domain_segment_app!(
            controller.app,
            controller.plot_grid,
            controller.selected_segment,
            controller.split_index;
            title_obs = controller.title_obs,
        )
        success || error("The selected split point is not valid.")
        update_partition_bindings!(controller; reset_index = true)
    end
end


function merge_boundary!(controller::QMLController, one_based_boundary)
    return guarded_action(controller, "Domain merge failed") do
        boundary = Int(one_based_boundary)
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
    end
end


function swap_boundary!(controller::QMLController, one_based_boundary)
    return guarded_action(controller, "Domain swap failed") do
        boundary = Int(one_based_boundary)
        selected_before_swap = controller.selected_segment
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
    end
end


function delete_selected_segment!(controller::QMLController)
    return guarded_action(controller, "Domain deletion failed") do
        deleted_segment = controller.selected_segment
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
    equation_images_by_model, equation_widths_by_model =
        render_equation_catalog(registry)
    labels = RD.model_labels(registry)
    isempty(labels) && error("No models found in MODEL_REGISTRY.")
    first_label = first(labels)
    first_key = RD.model_registry_key_for_label(registry, first_label)
    first_model = RD.get_model(registry, first_key)
    simulation = RD.create_simulation_state(
        first_model;
        N = N,
        dtmax = dtmax0,
        reltol = reltol,
        abstol = abstol,
        boundary_condition = boundary_condition0,
    )
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
    )
    app.plot_panel = RD.build_plot_panel!(plot_grid, app; title_obs = title)
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
    install_qml_renderfunction!(controller)
    register_qml_functions!(controller)
    properties = qml_property_map(
        controller,
        model_catalog_json(controller.registry),
        auto_close_ms = auto_close_ms,
    )
    QML.loadqml(QML_FILE; plot = figure, ui = properties)

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
