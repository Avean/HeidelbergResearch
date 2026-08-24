# Models/Mechanochemical/MechanochemicalGlobal.jl

# ============================================================
# Mechanochemical model with global normalization on the unit torus
# ============================================================
#
#     u_t = D u_xx - u + kappa exp(u(x)) / integral_T exp(u(y)) dy
#
# The periodic grid represents T = [0, 1), so its rectangle-rule weight is dx.
# ============================================================

RDModel(
    id = :mechanochemical_global,

    display_name = "Mechanochemical — global integral",

    variables = (:u,),

    parameters = (
        Du = 1e-3,
        kappa = 1.0,
    ),

    default_boundary_condition = :periodic,

    initial = function (U, x, p)
        U.u .= 0.01 .* randn(length(x))
        return nothing
    end,

    reaction = function (F, U, x, p, t)
        dx = x[2] - x[1]

        # A common exponential shift cancels between numerator and denominator
        # and prevents overflow for large positive values of u.
        shift = maximum(U.u)
        denominator = zero(eltype(U.u))

        for value in U.u
            denominator += exp(value - shift)
        end

        denominator *= dx

        for i in eachindex(U.u)
            F.u[i] = -U.u[i] +
                     p.kappa * exp(U.u[i] - shift) / denominator
        end

        return nothing
    end,

    diffusion = (
        u = :Du,
    ),

    latex_equations = (
        raw"\partial_t u = D\partial_{xx}u-u+\kappa\frac{e^{u(x)}}{\int_{\mathbb{T}}e^{u(y)}\,dy}",
    ),
)
