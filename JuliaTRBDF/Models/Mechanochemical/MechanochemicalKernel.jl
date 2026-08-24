# Models/Mechanochemical/MechanochemicalKernel.jl

# ============================================================
# Mechanochemical model with local normalization on the unit torus
# ============================================================
#
#     u_t = D u_xx - u
#           + kappa exp(u(x)) / integral_T K(x-y) exp(u(y)) dy
#
# K is the characteristic function of the periodic interval [x-0.1, x+0.1].
# Distances crossing 0 or 1 are wrapped modulo one.
# ============================================================

RDModel(
    id = :mechanochemical_kernel,

    display_name = "Mechanochemical — local kernel",

    variables = (:u,),

    parameters = (
        Du = 1e-5,
        kappa = 0.1, #0.4
    ),

    default_boundary_condition = :periodic,

    initial = function (U, x, p)
        U.u .= 0.01 .* randn(length(x))
        return nothing
    end,

    reaction = function (F, U, x, p, t)
        dx = x[2] - x[1]
        radius = 0.1 #0.3
        period = 1.0
        tolerance = 16.0 * eps(Float64)

        # Use a common shift in all exponentials. It cancels in every ratio.
        shift = maximum(U.u)
        shifted_exponential = similar(U.u)

        for i in eachindex(U.u)
            shifted_exponential[i] = exp(U.u[i] - shift)
        end

        for i in eachindex(U.u)
            denominator = zero(eltype(U.u))

            for j in eachindex(U.u)
                raw_distance = abs(x[i] - x[j])
                wrapped_distance = mod(raw_distance, period)
                torus_distance = min(
                    wrapped_distance,
                    period - wrapped_distance,
                )

                if torus_distance <= radius + tolerance
                    denominator += shifted_exponential[j]
                end
            end

            denominator *= dx
            F.u[i] = -U.u[i] +
                     p.kappa * shifted_exponential[i] / denominator
        end

        return nothing
    end,

    diffusion = (
        u = :Du,
    ),

    latex_equations = (
        raw"\partial_t u = D\partial_{xx}u-u+\kappa\frac{e^{u(x)}}{\int_{\mathbb{T}}K(x-y)e^{u(y)}\,dy}",
        raw"K(x-y)=\mathbf{1}_{\{d_{\mathbb{T}}(x,y)\leq 0.1\}}",
    ),
)
