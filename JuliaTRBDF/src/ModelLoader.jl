# src/ModelLoader.jl

# ============================================================
# Model registry
# ============================================================
#
# Model files are loaded once when the application module is loaded.
#
# Each file in models/ must evaluate to a ModelSpec object.
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
    # Return all Julia model files from the given directory.

    isdir(model_dir) ||
        error("Model directory does not exist: $model_dir")

    files = readdir(model_dir; join = true)
    files = filter(path -> endswith(path, ".jl"), files)

    return sort(files)
end


function model_file_labels(model_dir::AbstractString)
    # Return display labels for model files.

    files = model_files(model_dir)

    return basename.(files)
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
    #     Dict(filename => ModelSpec)

    registry = Dict{String, ModelSpec}()

    for path in model_files(model_dir)
        label = basename(path)
        registry[label] = load_model_from_file_at_startup(path)
    end

    isempty(registry) &&
        error("No model files found in directory: $model_dir")

    return registry
end


function model_labels(registry::Dict{String, ModelSpec})
    # Return sorted model labels.

    return sort(collect(keys(registry)))
end


function get_model(registry::Dict{String, ModelSpec}, label::String)
    # Get model by file label.

    haskey(registry, label) ||
        error("Unknown model label: $label")

    return registry[label]
end