# models/FisherKPP.jl

# ============================================================
# Fisher-KPP equation
# ============================================================
#
# Model:
#
#     u_t = D u_xx + u(1 - u)
#
# The diffusion part D u_xx is added automatically through:
#
#     diffusion = (
#         u = :D,
#     )
#
# Therefore the reaction function only defines:
#
#     F.u = u(1 - u)
#
# ============================================================

RDModel(
    id = :fisher_kpp,

    display_name = "Fisher-KPP",

    variables = (:u,),

    parameters = (
        D = 0.01,
    ),

    initial = function (U, x, p)
        Random.seed!(1)

        u = U.u

        @. u = 0.2 + 0.05 * cos(2π * x)
        u .+= 0.02 .* randn(length(x))

        return nothing
    end,

    reaction = function (F, U, x, p, t)
        
        @. F.u = U.u * (1.0 - U.u)

        return nothing
    end,

    diffusion = (
        u = :D,
    ),
)