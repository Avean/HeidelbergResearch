# models/SchnakenbergWithRho.jl

# ============================================================
# Schnakenberg system with active spatial rho profile
# ============================================================
#
# Model:
#
#     u_t = Du u_xx + ρ(x) (a - u + u^2 v)
#     v_t = Dv v_xx + b - u^2 v
#
# The active spatial profile ρ(x) is selected from the UI.
#
# In the reaction term, ρ is accessed as:
#
#     p.ρ
#
# The UI can switch between:
#
#     linear
#     gaussian
#     sinusoidal
#
# and this changes the RHS during simulation.
#
# ============================================================

RDModel(
    id = :schnakenberg_with_rho,

    display_name = "Schnakenberg with active rho(x)",

    variables = (:u, :v),

    parameters = (
        Du = 1e-4,
        Dv = 1e-2,

        a = 0.1,
        b = 0.9,

        ρ0 = 1.0,
        ρ1 = 0.5,
    ),

    initial = function (U, x, p)
        Random.seed!(5)

        u0 = p.a + p.b
        v0 = p.b / u0^2

        U.u .= u0 .+ 0.01 .* randn(length(x))
        U.v .= v0 .+ 0.01 .* randn(length(x))

        return nothing
    end,

    reaction = function (F, U, x, p, t)
        @. F.u = p.ρ * (p.a - U.u + U.u^2 * U.v)
        @. F.v = p.b - U.u^2 * U.v

        return nothing
    end,

    diffusion = (
        u = :Du,
        v = :Dv,
    ),

    spatial_profiles = (
        linear = (
            ρ = (x, p) -> (@. p.ρ0 + p.ρ1 * x),
        ),

        antilinear = (
            ρ = (x, p) -> (@. p.ρ0 + p.ρ1 * (1.0 - x)),
        ),

        sinusoidal = (
            ρ = (x, p) -> (@. p.ρ0 + p.ρ1 * sin(2π * x)),
        ),
    ),
)