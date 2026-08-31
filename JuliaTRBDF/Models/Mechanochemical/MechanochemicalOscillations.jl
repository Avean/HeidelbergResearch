# Models/Mechanochemical/MechanochemicalOscillations.jl

# ============================================================
# Mechanochemical model with a periodically reset linear ramp
# ============================================================
#
#     u_t = D u_xx - u
#           + kappa * ell(t) * exp(u(x)) / integral_T exp(u(y)) dy
#
# During every period of length T = 10, ell(t) grows linearly from zero
# towards 2.2 / kappa. Consequently, kappa * ell(t) grows from zero towards
# 2.2 and is reset to zero at the beginning of the next period.
# ============================================================

RDModel(
    id = :mechanochemical_oscillations,

    display_name = "Mechanochemical — oscillating ramp",

    variables = (:u,),

    parameters = (
        Du = 1e-3,
        kappa = 1.0,
        ramp_period = 10.0,
        ramp_threshold = 2.1,
    ),

    default_boundary_condition = :periodic,

    initial = function (U, x, p)
        U.u .= 0.01 .* randn(length(x))
        return nothing
    end,

    reaction = function (F, U, x, p, t)
        p.kappa != 0.0 || error("kappa must be nonzero in the oscillating-ramp model.")
        p.ramp_period > 0.0 || error("ramp_period must be positive.")

        dx = x[2] - x[1]
        phase = mod(t, p.ramp_period) / p.ramp_period
        ell = phase * (p.ramp_threshold / p.kappa)
        kappa_ell = p.kappa * ell

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
                     kappa_ell * exp(U.u[i] - shift) / denominator
        end

        return nothing
    end,

    diffusion = (
        u = :Du,
    ),

    latex_equations = (
        raw"\partial_t u=D\partial_{xx}u-u+\kappa\ell(t)\frac{e^{u(x)}}{\int_{\mathbb{T}}e^{u(y)}\,dy}",
        raw"\kappa\ell(t)=2.1\left( t/10 - [t/10] \right)",
    ),
)
