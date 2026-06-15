# src/ModelDSL.jl

# ============================================================
# Model DSL
# ============================================================
#
# User-facing layer for defining reaction-diffusion models.
#
# Supports:
#
#     - variables U.u, U.v, ...
#     - reaction fields F.u, F.v, ...
#     - parameters p.a, p.Du, ...
#     - spatial profiles p.ρ, p.source, ...
#     - active spatial profile sets switchable from the UI
#
# ============================================================


const ACTIVE_SPATIAL_PROFILE_SET_PARAM = :__active_spatial_profile_set_index


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
            nothing

        elseif coeff isa Function
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
# Spatial profile sets
# ============================================================

function _normalize_spatial_profiles(spatial_profiles)
    spatial_profiles isa NamedTuple ||
        error("spatial_profiles must be a NamedTuple.")

    return spatial_profiles
end


function _wrap_one_spatial_profile_set(
    profile_set::NamedTuple,
    param_names::Vector{Symbol},
)
    wrapped_profiles = Tuple{String, Function}[]

    for name in propertynames(profile_set)
        profile_fun = getproperty(profile_set, name)

        profile_fun isa Function ||
            error("Spatial profile $(name) must be a function.")

        wrapped_profile_fun = function (x, p_dict)
            p = _make_named_parameters(p_dict, param_names)
            return profile_fun(x, p)
        end

        push!(
            wrapped_profiles,
            (string(name), wrapped_profile_fun),
        )
    end

    return wrapped_profiles
end


function _is_flat_spatial_profile_tuple(spatial_profiles::NamedTuple)
    names = propertynames(spatial_profiles)

    isempty(names) &&
        return false

    return all(name -> getproperty(spatial_profiles, name) isa Function, names)
end


function _is_spatial_profile_set_tuple(spatial_profiles::NamedTuple)
    names = propertynames(spatial_profiles)

    isempty(names) &&
        return false

    return all(name -> getproperty(spatial_profiles, name) isa NamedTuple, names)
end


function _wrap_spatial_profile_sets(
    spatial_profiles::NamedTuple,
    param_names::Vector{Symbol},
)
    names = propertynames(spatial_profiles)

    if isempty(names)
        return Tuple{String, Vector{Tuple{String, Function}}}[]
    end

    # Backward-compatible flat syntax:
    #
    #     spatial_profiles = (
    #         ρ = (x, p) -> ...,
    #     )
    #
    # becomes one set named "default".
    if _is_flat_spatial_profile_tuple(spatial_profiles)
        wrapped_profiles = _wrap_one_spatial_profile_set(
            spatial_profiles,
            param_names,
        )

        return [
            ("default", wrapped_profiles),
        ]
    end

    # Grouped syntax:
    #
    #     spatial_profiles = (
    #         linear = (
    #             ρ = (x, p) -> ...,
    #         ),
    #
    #         gaussian = (
    #             ρ = (x, p) -> ...,
    #         ),
    #     )
    if _is_spatial_profile_set_tuple(spatial_profiles)
        wrapped_sets = Tuple{String, Vector{Tuple{String, Function}}}[]

        for set_name in names
            profile_set = getproperty(spatial_profiles, set_name)

            wrapped_profiles = _wrap_one_spatial_profile_set(
                profile_set,
                param_names,
            )

            push!(
                wrapped_sets,
                (string(set_name), wrapped_profiles),
            )
        end

        return wrapped_sets
    end

    error("""
    Invalid spatial_profiles format.

    Use either flat syntax:

        spatial_profiles = (
            ρ = (x, p) -> ...,
        )

    or grouped syntax:

        spatial_profiles = (
            linear = (
                ρ = (x, p) -> ...,
            ),
            gaussian = (
                ρ = (x, p) -> ...,
            ),
        )
    """)
end


function _active_spatial_profile_set_index(
    p_dict::Dict{Symbol, Float64},
    spatial_profile_sets,
)
    isempty(spatial_profile_sets) &&
        return 0

    if !haskey(p_dict, ACTIVE_SPATIAL_PROFILE_SET_PARAM)
        return 1
    end

    index = round(Int, p_dict[ACTIVE_SPATIAL_PROFILE_SET_PARAM])

    return clamp(index, 1, length(spatial_profile_sets))
end


function _evaluate_spatial_profile_for_parameter(
    x::AbstractVector,
    p_dict::Dict{Symbol, Float64},
    profile_name::String,
    profile_fun::Function,
)
    raw = profile_fun(x, p_dict)

    y = if raw isa Number
        fill(Float64(raw), length(x))
    else
        Float64.(collect(raw))
    end

    length(y) == length(x) ||
        error("Spatial profile $(profile_name) has wrong length.")

    return y
end


function _make_named_parameters_with_active_spatial_profiles(
    p_dict::Dict{Symbol, Float64},
    param_names::Vector{Symbol},
    spatial_profile_sets,
    x::AbstractVector,
)
    p_base = _make_named_parameters(p_dict, param_names)

    isempty(spatial_profile_sets) &&
        return p_base

    set_index = _active_spatial_profile_set_index(
        p_dict,
        spatial_profile_sets,
    )

    _, profiles = spatial_profile_sets[set_index]

    profile_symbols = Symbol[]
    profile_values = Any[]

    for (profile_name, profile_fun) in profiles
        profile_symbol = Symbol(profile_name)

        profile_symbol in param_names &&
            error("Spatial profile $(profile_name) conflicts with an existing parameter name.")

        push!(profile_symbols, profile_symbol)

        y = _evaluate_spatial_profile_for_parameter(
            x,
            p_dict,
            profile_name,
            profile_fun,
        )

        push!(profile_values, y)
    end

    p_profiles = NamedTuple{Tuple(profile_symbols)}(Tuple(profile_values))

    return merge(p_base, p_profiles)
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

    wrapped_spatial_profile_sets = _wrap_spatial_profile_sets(
        spatial_profiles,
        param_names,
    )

    if !isempty(wrapped_spatial_profile_sets)
        default_params[ACTIVE_SPATIAL_PROFILE_SET_PARAM] = 1.0
    end

    function initialize_wrapped!(
        Umat::AbstractMatrix,
        x::AbstractVector,
        p_dict::Dict{Symbol, Float64},
    )
        size(Umat, 2) == length(vars) ||
            error("Initial condition matrix has wrong number of variables.")

        U = VariableViews(Umat, var_index)

        p = _make_named_parameters_with_active_spatial_profiles(
            p_dict,
            param_names,
            wrapped_spatial_profile_sets,
            x,
        )

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

        p = _make_named_parameters_with_active_spatial_profiles(
            p_dict,
            param_names,
            wrapped_spatial_profile_sets,
            x,
        )

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
        spatial_profile_sets = wrapped_spatial_profile_sets,
    )
end