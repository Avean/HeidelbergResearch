# src/ModelLoader.jl

# ============================================================
# Model registry
# ============================================================
#
# Model files are loaded once when the application module is loaded.
#
# Every model lives in Models/<family>/<model>/model.jl.  The first folder is
# the family used by the first UI selector and the second folder is the
# concrete model.  TOML files beside model.jl are Series presets for it.
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
    # Return model entry points from the supported three-level layout.

    isdir(model_dir) ||
        error("Model directory does not exist: $model_dir")

    files = String[]

    for family in readdir(model_dir; join = true)
        isdir(family) || continue
        for model_directory in readdir(family; join = true)
            isdir(model_directory) || continue
            path = joinpath(model_directory, "model.jl")
            isfile(path) && push!(files, path)
        end
    end

    return sort(files)
end


function model_registry_key(
    model_dir::AbstractString,
    path::AbstractString,
)
    relative = replace(relpath(dirname(path), model_dir), '\\' => '/')
    parts = split(relative, '/')
    length(parts) == 2 || error("Model path must be Models/<family>/<model>/model.jl: $path")
    return relative
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
    return words_from_identifier(basename(String(key)))
end


function model_menu_label(key::AbstractString, model::ModelSpec)
    # Families already have their own selector, so avoid repeating their name
    # in every entry of the concrete-model selector.
    label = model.display_name
    parts = split(label, '—'; limit = 2)
    length(parts) == 2 && (label = strip(parts[2]))

    # Some legacy display names use a plain prefix rather than an em dash,
    # e.g. "MathBio basic".  Strip that family prefix only when it is a whole
    # first word, so model names themselves stay intact.
    family_compact = lowercase(replace(model_family_name(key), r"[\s_-]" => ""))
    words = split(label)
    if !isempty(words) &&
       lowercase(replace(first(words), r"[\s_-]" => "")) == family_compact
        label = join(words[2:end], " ")
    end

    return isempty(label) ? model.display_name : label
end


function model_menu_catalog(registry::Dict{String, ModelSpec})
    families = Dict{String, Vector{NamedTuple}}()

    for key in sort(collect(keys(registry)))
        family = model_family_name(key)
        entry = (
            key = key,
            label = model_menu_label(key, registry[key]),
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
