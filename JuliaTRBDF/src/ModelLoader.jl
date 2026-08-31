# src/ModelLoader.jl

# ============================================================
# Model registry
# ============================================================
#
# Model files are loaded once when the application module is loaded.
#
# Each Julia file anywhere below models/ must evaluate to a ModelSpec object.
# First-level folders are used as model families in hierarchical UI menus.
#
# Example:
#
#     ModelSpec(
#         id = :fisher_kpp,
#         ...
#     )
#
# The last expression in the file must be the ModelSpec.
#
# ============================================================


function model_files(model_dir::AbstractString)
    # Return all Julia model files recursively from the model directory.

    isdir(model_dir) ||
        error("Model directory does not exist: $model_dir")

    files = String[]

    for (root, _, names) in walkdir(model_dir)
        for name in names
            endswith(lowercase(name), ".jl") || continue
            push!(files, joinpath(root, name))
        end
    end

    return sort(files)
end


function model_registry_key(
    model_dir::AbstractString,
    path::AbstractString,
)
    return replace(relpath(path, model_dir), '\\' => '/')
end


function model_file_labels(model_dir::AbstractString)
    # Return display labels for model files.

    files = model_files(model_dir)

    return [model_registry_key(model_dir, path) for path in files]
end


function load_model_from_file_at_startup(path::AbstractString)::ModelSpec
    # Load one model file.
    #
    # Important:
    # This function is intended to be called while the main application
    # module is being loaded, not during the simulation loop.

    isfile(path) ||
        error("Model file does not exist: $path")

    model = Base.include(@__MODULE__, path)

    model isa ModelSpec ||
        error("""
        Model file must evaluate to a ModelSpec object: $path

        The file should end with something like:

            ModelSpec(
                id = :my_model,
                ...
            )
        """)

    validate_model(model)

    return model
end


function load_model_registry(model_dir::AbstractString)
    # Load all models from the models/ directory.
    #
    # Returns:
    #
    #     Dict(relative/path/filename => ModelSpec)

    registry = Dict{String, ModelSpec}()

    for path in model_files(model_dir)
        key = model_registry_key(model_dir, path)
        registry[key] = load_model_from_file_at_startup(path)
    end

    isempty(registry) &&
        error("No model files found in directory: $model_dir")

    return registry
end


function model_labels(registry::Dict{String, ModelSpec})
    # Preserve the old flat selector by showing filenames when they are unique.

    labels = basename.(collect(keys(registry)))

    if length(unique(labels)) != length(labels)
        return sort(collect(keys(registry)))
    end

    return sort(labels)
end


function get_model(registry::Dict{String, ModelSpec}, label::String)
    # Accept either the full registry key or an unambiguous filename.

    if haskey(registry, label)
        return registry[label]
    end

    matching_keys = filter(key -> basename(key) == label, keys(registry))

    isempty(matching_keys) && error("Unknown model label: $label")
    length(matching_keys) == 1 ||
        error("Ambiguous model filename; use its folder-qualified key: $label")

    return registry[only(matching_keys)]
end


function model_registry_key_for_label(
    registry::Dict{String, ModelSpec},
    label::String,
)
    haskey(registry, label) && return label

    matching_keys = filter(key -> basename(key) == label, keys(registry))
    length(matching_keys) == 1 ||
        error("Unknown or ambiguous model label: $label")

    return only(matching_keys)
end


function words_from_identifier(value::AbstractString)
    words = replace(String(value), '_' => ' ', '-' => ' ')
    words = replace(words, r"(?<=[a-z0-9])(?=[A-Z])" => " ")
    words = replace(words, r"(?<=[A-Za-z])(?=[0-9])" => " ")

    return strip(words)
end


function model_family_name(key::AbstractString)
    directory = dirname(String(key))
    directory == "." && return "Other"

    return words_from_identifier(first(split(replace(directory, '\\' => '/'), '/')))
end


function model_variant_name(key::AbstractString)
    stem = splitext(basename(String(key)))[1]
    directory_stem = basename(dirname(String(key)))

    if directory_stem != "." && startswith(stem, directory_stem)
        remainder = stem[(length(directory_stem) + 1):end]
        isempty(remainder) && return "Basic"

        return words_from_identifier(remainder)
    end

    return words_from_identifier(stem)
end


function model_menu_catalog(registry::Dict{String, ModelSpec})
    families = Dict{String, Vector{NamedTuple}}()

    for key in sort(collect(keys(registry)))
        family = model_family_name(key)
        entry = (
            key = key,
            label = model_variant_name(key),
        )
        push!(get!(families, family, NamedTuple[]), entry)
    end

    return [
        (
            family = family,
            models = sort(families[family]; by = entry -> entry.label),
        )
        for family in sort(collect(keys(families)))
    ]
end
