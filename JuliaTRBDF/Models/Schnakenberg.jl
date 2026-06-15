# models/Schnakenberg.jl

# ============================================================
# Schnakenberg reaction-diffusion system
# ============================================================
#
# Model:
#
#     u_t = Du u_xx + a - u + u^2 v
#     v_t = Dv v_xx + b - u^2 v
#
# Homogeneous steady state:
#
#     u_* = a + b
#     v_* = b / (a + b)^2
#
# The diffusion terms are added automatically through:
#
#     diffusion = (
#         u = :Du,
#         v = :Dv,
#     )
#
# Therefore the reaction function only defines:
#
#     F.u = a - u + u^2 v
#     F.v = b - u^2 v
#
# ============================================================

RDModel(
    id = :schnakenberg,

    display_name = "Schnakenberg system",

    variables = (:u, :v),

    parameters = (
        Du = 1e-4,
        Dv = 1e-2,
        a  = 0.1,
        b  = 0.9,
    ),

    initial = function (U, x, p)
        Random.seed!(2)
        
        u0 = p.a + p.b
        v0 = p.b / u0^2

        U.u .= u0 .+ 0.01 .* randn(length(x))
        U.v .= v0 .+ 0.01 .* randn(length(x))

        return nothing
    end,

    reaction = function (F, U, x, p, t)
        
        @. F.u = p.a - U.u + U.u^2 * U.v
        @. F.v = p.b - U.u^2 * U.v

        return nothing
    end,

    diffusion = (
        u = :Du,
        v = :Dv,
    ),
)