# Models/Mechanochemical/MechanochemicalKernel2.jl

# ============================================================
# Mechanochemical model with globally normalized local interaction
# ============================================================
#
# Let
#
#     A(x) = integral_T K(x,y) exp(u(y)) dy,
#     Z    = integral_T exp(u(x)) A(x) dx.
#
# The equation is
#
#     u_t = D u_xx - u + kappa exp(u(x)) A(x) / Z.
#
# K grows linearly with periodic physical distance d:
#
#     K(d) = 1 + 4d  for d <= 1/4,
#     K(d) = 0       for d > 1/4.
#
# The physical cutoff 1/4 remains fixed when the domain is rescaled.
# ============================================================

RDModel(
    id = :mechanochemical_kernel_2,

    display_name = "Mechanochemical — normalized kernel",

    variables = (:u,),

    parameters = (
        Du = 1e-4,
        kappa = 1.0,
        domain_scale = 1.0,
    ),

    default_boundary_condition = :periodic,

    initial = function (U, x, p)
        U.u .= 0.01 .* randn(length(x))
        return nothing
    end,

    reaction = function (F, U, x, p, t)
        base_dx = x[2] - x[1]
        physical_dx = p.domain_scale * base_dx
        physical_cutoff = 1.0 / 4.0
        base_period = 1.0
        tolerance = 16.0 * eps(Float64)

        # A common exponential shift cancels from both the numerator and the
        # double-integral normalization, preventing overflow for large u.
        shift = maximum(U.u)
        shifted_exponential = similar(U.u)
        local_integral = similar(U.u)

        for i in eachindex(U.u)
            shifted_exponential[i] = exp(U.u[i] - shift)
        end

        for i in eachindex(U.u)
            integral_value = zero(eltype(U.u))

            for j in eachindex(U.u)
                raw_distance = abs(x[i] - x[j])
                wrapped_distance = mod(raw_distance, base_period)
                base_torus_distance = min(
                    wrapped_distance,
                    base_period - wrapped_distance,
                )
                physical_distance = p.domain_scale * base_torus_distance

                if physical_distance <= physical_cutoff + tolerance
                    kernel_value = 1.0 + 4.0 * physical_distance
                    integral_value += kernel_value * shifted_exponential[j]
                end
            end

            local_integral[i] = physical_dx * integral_value
        end

        normalization = zero(eltype(U.u))

        for i in eachindex(U.u)
            normalization += shifted_exponential[i] * local_integral[i]
        end

        normalization *= physical_dx

        for i in eachindex(U.u)
            nonlocal_term = shifted_exponential[i] *
                            local_integral[i] / normalization
            F.u[i] = -U.u[i] + p.kappa * nonlocal_term
        end

        return nothing
    end,

    diffusion = (
        u = :Du,
    ),

    latex_equations = (
        raw"\partial_t u = D\partial_{xx}u-u+\kappa\frac{e^{u(x)}\int_{\mathbb{T}}K(x,y)e^{u(y)}\,dy}{\int_{\mathbb{T}}\int_{\mathbb{T}}e^{u(x)}e^{u(y)}K(x,y)\,dx\,dy}",
        raw"K(x,y)=\left(1+4d_{\mathbb{T}}(x,y)\right)\mathbf{1}_{\left\{d_{\mathbb{T}}(x,y)\leq\frac{1}{4}\right\}}",
    ),
)
