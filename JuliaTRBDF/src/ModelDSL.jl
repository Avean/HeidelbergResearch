# src/ModelDSL.jl

# ============================================================
# User-friendly model definition layer
# ============================================================
#
# This file defines RDModel(...), a more convenient interface for
# writing models.
#
# The application still uses ModelSpec internally.
# RDModel(...) simply translates a user-friendly model definition
# into a standard ModelSpec.
#
# User model files can write:
#
#     RDModel(
#         id = :fisher_kpp,
#         display_name = "Fisher-KPP",
#         variables = (:u,),
#         parameters = (
#             D = 0.01,
#         ),
#         initial = function (U, x, p)
#             @. U.u = 0.2 + 0.05 * cos(2π * x)
#         end,
#         reaction = function (F, U, x, p, t)
#             @. F.u = U.u * (1.0 - U.u)
#         end,
#         diffusion = (
#             u = :D,
#         ),
#     )
#
# The user does not need to write:
#
#     U[:, 1]
#     dU[:, 1]
#     p[:D]
#     mul!(du, Lap, u)
#
# ============================================================


# ============================================================
# Named variable views
# ============================================================

struct VariableViews{M <: AbstractMatrix}
    data::M
    var_index::Dict{Symbol, Int}
end


function Base.getproperty(V::VariableViews, name::Symbol)
    # Allow access like:
    #
    #     U.u
    #     U.v
    #     F.rho

    if name === :data || name === :var_index
        return getfield(V, name)
    end

    var_index = getfield(V, :var_index)

    if haskey(var_index, name)
        j = var_index[name]
        return view(getfield(V, :data), :, j)
    end

    error("Unknown model variable: $name")
end


function Base.propertynames(V::VariableViews; private::Bool = false)
    names = collect(keys(getfield(V, :var_index)))

    if private
        return Tuple(vcat(names, [:data, :var_index]))
    else
        return Tuple(names)
    end
end


# ============================================================
# Parameter helpers
# ============================================================

function _make_parameter_dict(parameters::NamedTuple)
    # Convert a user-friendly NamedTuple of parameters into the internal
    # Dict{Symbol, Float64} representation expected by ModelSpec.

    params = Dict{Symbol, Float64}()

    for name in keys(parameters)
        value = getproperty(parameters, name)

        value isa Real ||
            error("Parameter $name must be a real number.")

        params[name] = Float64(value)
    end

    return params
end


function _make_named_parameters(
    params::Dict{Symbol, Float64},
    param_names::Vector{Symbol},
)
    # Convert internal Dict parameters back into a NamedTuple.
    #
    # This lets user model code write:
    #
    #     p.D
    #     p.κ
    #
    # instead of:
    #
    #     p[:D]
    #     p[:κ]

    return NamedTuple{Tuple(param_names)}(
        Tuple(params[name] for name in param_names)
    )
end


# ============================================================
# Variable and diffusion helpers
# ============================================================

function _normalize_variables(variables)
    # Accept either:
    #
    #     :u
    #     (:u,)
    #     (:u, :v)
    #     [:u, :v]

    if variables isa Symbol
        vars = [variables]

    elseif variables isa Tuple
        vars = collect(variables)

    elseif variables isa AbstractVector
        vars = collect(variables)

    else
        error("variables must be a Symbol, Tuple, or Vector of Symbols.")
    end

    all(v -> v isa Symbol, vars) ||
        error("All variable names must be Symbols.")

    length(vars) >= 1 ||
        error("A model must have at least one variable.")

    length(unique(vars)) == length(vars) ||
        error("Variable names must be unique.")

    return Symbol.(vars)
end


function _normalize_diffusion(diffusion)
    # Diffusion should usually be a NamedTuple, for example:
    #
    #     diffusion = (
    #         u = :D,
    #         v = :Dv,
    #     )
    #
    # A missing variable means zero diffusion.

    diffusion === nothing && return (;)

    diffusion isa NamedTuple ||
        error("diffusion must be a NamedTuple or nothing.")

    return diffusion
end


function _diffusion_coefficient(
    diffusion::NamedTuple,
    variable::Symbol,
    p,
)
    # Return the diffusion coefficient for one variable.
    #
    # Allowed forms:
    #
    #     diffusion = (
    #         u = :D,
    #     )
    #
    # means use parameter p.D.
    #
    #     diffusion = (
    #         u = 0.01,
    #     )
    #
    # means use constant 0.01.
    #
    # If the variable is missing from diffusion, the coefficient is 0.0.

    if !(variable in keys(diffusion))
        return 0.0
    end

    spec = getproperty(diffusion, variable)

    if spec isa Symbol
        return Float64(getproperty(p, spec))

    elseif spec isa Real
        return Float64(spec)

    else
        error("Invalid diffusion specification for variable $variable.")
    end
end


# ============================================================
# Main user-facing constructor
# ============================================================

function RDModel(;
    id::Symbol,
    display_name::String,
    variables,
    parameters::NamedTuple = (;),
    initial::Function,
    reaction::Function,
    diffusion = (;),
)
    # Build a standard internal ModelSpec from a user-friendly model
    # definition.
    #
    # User-facing convention:
    #
    #     U.u  = current value of variable u
    #     F.u  = reaction part of u_t
    #     p.D  = parameter D
    #
    # The diffusion part is added automatically from the diffusion field.

    vars = _normalize_variables(variables)
    diffusion = _normalize_diffusion(diffusion)

    param_names = collect(keys(parameters))
    default_params = _make_parameter_dict(parameters)

    var_index = Dict{Symbol, Int}()

    for (j, var) in enumerate(vars)
        var_index[var] = j
    end

    for var in keys(diffusion)
        var in vars ||
            error("Diffusion specified for unknown variable: $var")
    end

    initialize_wrapped! = function (Umat, x, p_dict)
        U = VariableViews(Umat, var_index)
        p = _make_named_parameters(p_dict, param_names)

        initial(U, x, p)

        return nothing
    end

    rhs_wrapped! = function (dUmat, Umat, Lap, x, p_dict, t)
        U = VariableViews(Umat, var_index)
        F = VariableViews(dUmat, var_index)
        p = _make_named_parameters(p_dict, param_names)

        # First compute the reaction part.
        #
        # User code writes:
        #
        #     F.u = reaction term
        #
        # so we start from zero.
        fill!(dUmat, zero(eltype(dUmat)))

        reaction(F, U, x, p, t)

        # Then add diffusion automatically:
        #
        #     u_t = D_u Lap u + F_u

        tmp = similar(view(Umat, :, 1))

        for var in vars
            j = var_index[var]

            D = _diffusion_coefficient(diffusion, var, p)

            if D != 0.0
                u  = view(Umat, :, j)
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
        varnames = String.(vars),
        default_params = default_params,
        initialize! = initialize_wrapped!,
        rhs! = rhs_wrapped!,
    )
end