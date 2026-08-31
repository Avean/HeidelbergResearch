# Models/Mechanochemical/MechanochemicalKernel3.jl

# ============================================================
# Mechanochemical model with a normalized ring kernel
# ============================================================
#
# Let
#
#     K0(x,y) = 10  for 1/4 <= d_T(x,y) <= 1/2 - eps,
#     K0(x,y) = 1   otherwise,
#
# and normalize each row so that
#
#     integral_T K(x,y) dy = 1.
#
# With
#
#     A(x) = integral_T K(x,y) exp(u(y)) dy,
#     Z    = integral_T exp(u(x)) A(x) dx,
#
# the equation is
#
#     u_t = D u_xx - u + kappa exp(u(x)) A(x) / Z.
#
# Both ring radii are physical and remain fixed when the domain is rescaled.
# ============================================================

RDModel(
    id = :mechanochemical_kernel_3,

    display_name = "Mechanochemical — normalized ring kernel",

    variables = (:u,),

    parameters = (
        Du = 1e-3,
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
        inner_radius = 1.0 / 4.0
        outer_radius = 1.0 / 2.0 - eps(Float64)
        base_period = 1.0
        tolerance = 16.0 * eps(Float64)

        # The common shift cancels from the final normalized expression and
        # prevents overflow for large positive values of u.
        shift = maximum(U.u)
        shifted_exponential = similar(U.u)
        local_integral = similar(U.u)

        for i in eachindex(U.u)
            shifted_exponential[i] = exp(U.u[i] - shift)
        end

        for i in eachindex(U.u)
            kernel_mass = zero(eltype(U.u))
            weighted_value = zero(eltype(U.u))

            for j in eachindex(U.u)
                raw_distance = abs(x[i] - x[j])
                wrapped_distance = mod(raw_distance, base_period)
                base_torus_distance = min(
                    wrapped_distance,
                    base_period - wrapped_distance,
                )
                physical_distance = p.domain_scale * base_torus_distance
                in_ring =
                    physical_distance + tolerance >= inner_radius &&
                    physical_distance <= outer_radius
                raw_kernel = in_ring ? 10.0 : 1.0

                kernel_mass += raw_kernel
                weighted_value += raw_kernel * shifted_exponential[j]
            end

            # Dividing the weighted integral by the raw kernel mass is
            # equivalent to first normalizing K so that integral K(x,y)dy = 1.
            local_integral[i] = weighted_value / kernel_mass
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
        raw"K_0(x,y)=1+9\mathbf{1}_{\left\{\frac{1}{4}\leq d_{\mathbb{T}}(x,y)\leq\frac{1}{2}-\varepsilon\right\}}",
        raw"K(x,y)=\frac{K_0(x,y)}{\int_{\mathbb{T}}K_0(x,z)\,dz}",
    ),
)
