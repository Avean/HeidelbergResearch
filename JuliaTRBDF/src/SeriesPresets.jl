# src/SeriesPresets.jl

# ============================================================
# File-backed presets for Series perturbations
# ============================================================

struct SeriesPresetPerturbation
    variable_name::String
    position::Symbol
    width_min_fraction::Float64
    width_max_fraction::Float64
    height_min::Float64
    height_max::Float64
end


struct SeriesPresetDefinition
    key::String
    name::String
    model_key::String
    initial_values::Dict{String, Float64}
    perturbations::Vector{SeriesPresetPerturbation}
end


function _series_preset_string(table, field::String, path::AbstractString)
    value = get(table, field, nothing)
    value isa AbstractString || error("Preset $path must define string field '$field'.")
    return String(value)
end


function _series_preset_number(table, field::String, path::AbstractString)
    value = get(table, field, nothing)
    value isa Real || error("Preset $path must define numeric field '$field'.")
    number = Float64(value)
    isfinite(number) || error("Preset $path field '$field' must be finite.")
    return number
end


function _series_preset_position(value::AbstractString, path::AbstractString)
    position = Symbol(lowercase(strip(value)))
    position in (:left, :center, :right) ||
        error("Preset $path position must be left, center or right.")
    return position
end


function load_series_preset(path::AbstractString, model_key::AbstractString)
    data = TOML.parsefile(path)
    name = _series_preset_string(data, "name", path)
    raw_initial_values = get(data, "initial_values", Dict{String, Any}())
    raw_initial_values isa AbstractDict ||
        error("Preset $path initial_values must be a table.")
    initial_values = Dict{String, Float64}()
    for (variable, value) in raw_initial_values
        value isa Real || error("Preset $path initial value for $variable must be numeric.")
        number = Float64(value)
        isfinite(number) || error("Preset $path initial value for $variable must be finite.")
        initial_values[String(variable)] = number
    end
    raw_perturbations = get(data, "perturbations", nothing)
    raw_perturbations isa AbstractVector ||
        error("Preset $path must define at least one [[perturbations]] entry.")
    isempty(raw_perturbations) && error("Preset $path has no perturbations.")

    perturbations = SeriesPresetPerturbation[]
    for entry in raw_perturbations
        entry isa AbstractDict || error("Preset $path contains an invalid perturbation.")
        variable_name = _series_preset_string(entry, "variable", path)
        position = _series_preset_position(_series_preset_string(entry, "position", path), path)
        width_min = _series_preset_number(entry, "width_min_fraction", path)
        width_max = _series_preset_number(entry, "width_max_fraction", path)
        height_min = _series_preset_number(entry, "height_min", path)
        height_max = _series_preset_number(entry, "height_max", path)
        0.0 < width_min <= width_max ||
            error("Preset $path has invalid width bounds.")
        height_min <= height_max || error("Preset $path has invalid height bounds.")
        push!(
            perturbations,
            SeriesPresetPerturbation(
                variable_name,
                position,
                width_min,
                width_max,
                height_min,
                height_max,
            ),
        )
    end

    key = basename(path)
    return SeriesPresetDefinition(key, name, String(model_key), initial_values, perturbations)
end


function load_series_presets(root::AbstractString)
    isdir(root) || return SeriesPresetDefinition[]
    presets = SeriesPresetDefinition[]

    for model_path in RD.model_files(root)
        directory = dirname(model_path)
        model_key = RD.model_registry_key(root, model_path)
        for name in sort(readdir(directory))
            endswith(lowercase(name), ".toml") || continue
            path = joinpath(directory, name)
            try
                push!(presets, load_series_preset(path, model_key))
            catch error
                @warn "Ignoring invalid Series preset" path exception = (error, catch_backtrace())
            end
        end
    end

    sort!(presets; by = preset -> (preset.model_key, preset.name, preset.key))
    return presets
end


function series_presets_json(presets::Vector{SeriesPresetDefinition}, model_key::AbstractString)
    entries = ["{\"key\":\"none\",\"name\":\"None\"}"]

    for preset in presets
        preset.model_key == model_key || continue
        push!(
            entries,
            "{\"key\":" * json_string(preset.key) *
            ",\"name\":" * json_string(preset.name) * "}",
        )
    end

    return "[" * join(entries, ",") * "]"
end
