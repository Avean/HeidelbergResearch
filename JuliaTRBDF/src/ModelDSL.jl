# src/ModelDSL.jl

# ============================================================
# Model DSL
# ============================================================
#
# User-facing layer for defining reaction-diffusion models.
#
# A model file can end with:
#
#     RDModel(
#         id = :my_model,
#         display_name = "My model",
#         variables = (:u, :v),
#         parameters = (...),
#         initial = function (U, x, p)
#             ...
#         end,
#         reaction = function (F, U, x, p, t)
#             ...
#         end,
#         diffusion = (
#             u = :Du,
#             v = :Dv,
#         ),
#         spatial_profiles = (
#             ρ = (x, p) -> @. p.ρ0 + p.ρ1 * x,
#         ),
#     )
#
# The user writes:
#
#     U.u, U.v       current variables
#     F.u, F.v       reaction right-hand sides
#     p.Du, p.a      parameters
#     x              spatial grid
#
# Diffusion is added automatically by the DSL.
#
# ============================================================


# ============================================================
# Variable views
# ============================================================

struct VariableViews{M <: AbstractMatrix}
    data::M
    var_index::Dict{Symbol, Int}
end


function Base.getproperty(V::VariableViews, name::Symbol)
    if name === :data || name === :var_index
        return getfield(V, name)
    end

    if haskey(V.var_index, name)
        return view(V.data, :, V.var_index[name])
    end

    error("Unknown variable: $(name)")
end


function Base.propertynames(V::VariableViews; private::Bool = false)
    names = collect(keys(V.var_index))

    if private
        append!(names, [:data, :var_index])
    end

    return Tuple(names)
end


# ============================================================
# Parameters
# ============================================================

function _make_parameter_dict(parameters::NamedTuple)
    params = Dict{Symbol, Float64}()

    for name in propertynames(parameters)
        value = getproperty(parameters, name)

        value isa Real ||
            error("Parameter $(name) must be a real number.")

        params[name] = Float64(value)
    end

    return params
end


function _make_named_parameters(
    params::Dict{Symbol, Float64},
    param_names::Vector{Symbol},
)
    values = map(name -> params[name], param_names)
    return NamedTuple{Tuple(param_names)}(Tuple(values))
end


# ============================================================
# Variables
# ============================================================

function _normalize_variables(variables)
    vars = collect(Symbol.(variables))

    isempty(vars) &&
        error("Model must contain at least one variable.")

    length(unique(vars)) == length(vars) ||
        error("Variable names must be unique.")

    return vars
end


function _make_variable_index(vars::Vector{Symbol})
    var_index = Dict{Symbol, Int}()

    for (j, var) in enumerate(vars)
        var_index[var] = j
    end

    return var_index
end


# ============================================================
# Diffusion
# ============================================================

function _normalize_diffusion(diffusion)
    diffusion isa NamedTuple ||
        error("diffusion must be a NamedTuple, for example (u = :Du, v = :Dv).")

    return diffusion
end


function _validate_diffusion(
    diffusion::NamedTuple,
    vars::Vector{Symbol},
    param_names::Vector{Symbol},
)
    for var in propertynames(diffusion)
        var in vars ||
            error("Diffusion coefficient specified for unknown variable $(var).")

        coeff = getproperty(diffusion, var)

        if coeff isa Symbol
            coeff in param_names ||
                error("Unknown diffusion parameter $(coeff) for variable $(var).")

        elseif coeff isa Real
            # Constant numeric diffusion coefficient.
            nothing

        elseif coeff isa Function
            # Function of parameters, e.g.
            #
            #     diffusion = (
            #         u = p -> p.Du,
            #     )
            nothing

        else
            error("""
            Invalid diffusion coefficient for variable $(var).

            Allowed forms are:

                u = :Du
                u = 0.01
                u = p -> p.Du
            """)
        end
    end

    return nothing
end


function _diffusion_coefficient(
    diffusion::NamedTuple,
    var::Symbol,
    p,
)
    if !(var in propertynames(diffusion))
        return 0.0
    end

    coeff = getproperty(diffusion, var)

    if coeff isa Symbol
        return getproperty(p, coeff)

    elseif coeff isa Real
        return Float64(coeff)

    elseif coeff isa Function
        return coeff(p)

    else
        error("Invalid diffusion coefficient for variable $(var).")
    end
end


# ============================================================
# Spatial profiles
# ============================================================

function _normalize_spatial_profiles(spatial_profiles)
    spatial_profiles isa NamedTuple ||
        error("spatial_profiles must be a NamedTuple.")

    return spatial_profiles
end


function _wrap_spatial_profiles(
    spatial_profiles::NamedTuple,
    param_names::Vector{Symbol},
)
    wrapped_spatial_profiles = Tuple{String, Function}[]

    for name in propertynames(spatial_profiles)
        profile_fun = getproperty(spatial_profiles, name)

        profile_fun isa Function ||
            error("Spatial profile $(name) must be a function.")

        wrapped_profile_fun = function (x, p_dict)
            p = _make_named_parameters(p_dict, param_names)
            return profile_fun(x, p)
        end

        push!(
            wrapped_spatial_profiles,
            (string(name), wrapped_profile_fun),
        )
    end

    return wrapped_spatial_profiles
end


# ============================================================
# Main model constructor
# ============================================================

function RDModel(;
    id::Symbol,
    display_name::String,
    variables,
    parameters,
    initial::Function,
    reaction::Function,
    diffusion = NamedTuple(),
    spatial_profiles = NamedTuple(),
)
    vars = _normalize_variables(variables)
    var_index = _make_variable_index(vars)

    parameters isa NamedTuple ||
        error("parameters must be a NamedTuple.")

    default_params = _make_parameter_dict(parameters)
    param_names = collect(Symbol.(propertynames(parameters)))

    diffusion = _normalize_diffusion(diffusion)
    _validate_diffusion(diffusion, vars, param_names)

    spatial_profiles = _normalize_spatial_profiles(spatial_profiles)

    wrapped_spatial_profiles = _wrap_spatial_profiles(
        spatial_profiles,
        param_names,
    )

    function initialize_wrapped!(
        Umat::AbstractMatrix,
        x::AbstractVector,
        p_dict::Dict{Symbol, Float64},
    )
        size(Umat, 2) == length(vars) ||
            error("Initial condition matrix has wrong number of variables.")

        U = VariableViews(Umat, var_index)
        p = _make_named_parameters(p_dict, param_names)

        initial(U, x, p)

        return nothing
    end

    function rhs_wrapped!(
        dUmat::AbstractMatrix,
        Umat::AbstractMatrix,
        Lap,
        x::AbstractVector,
        p_dict::Dict{Symbol, Float64},
        t,
    )
        size(Umat, 2) == length(vars) ||
            error("State matrix has wrong number of variables.")

        size(dUmat, 2) == length(vars) ||
            error("RHS matrix has wrong number of variables.")

        U = VariableViews(Umat, var_index)
        F = VariableViews(dUmat, var_index)
        p = _make_named_parameters(p_dict, param_names)

        fill!(dUmat, zero(eltype(dUmat)))

        reaction(F, U, x, p, t)

        tmp = similar(view(Umat, :, 1))

        for var in vars
            j = var_index[var]
            D = _diffusion_coefficient(diffusion, var, p)

            if D != 0.0
                u = view(Umat, :, j)
                du = view(dUmat, :, j)

                mul!(tmp, Lap, u)

                @. du += D * tmp
            end
        end

        return nothing
    end

    return ModelSpec(
        id = id,
        display_name = display_name,
        nvars = length(vars),
        varnames = string.(vars),
        default_params = default_params,
        initialize! = initialize_wrapped!,
        rhs! = rhs_wrapped!,
        spatial_profiles = wrapped_spatial_profiles,
    )
end