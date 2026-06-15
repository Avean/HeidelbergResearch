# models/FisherKPP.jl

# ============================================================
# Fisher-KPP equation
# ============================================================
#
# Model:
#
#     u_t = D u_xx + u(1 - u)
#
# ============================================================

ModelSpec(
    id = :fisher_kpp,

    display_name = "Fisher-KPP",

    nvars = 1,

    varnames = ["u"],

    default_params = Dict(
        :D => 0.01,
    ),

    initialize! = function (U, x, p)
        Random.seed!(1)

        N = length(x)

        U[:, 1] .= @. 0.2 + 0.05 * cos(2π * x)
        U[:, 1] .+= 0.02 .* randn(N)

        return nothing
    end,

    rhs! = function (dU, U, Lap, x, p, t)
        u  = @view U[:, 1]
        du = @view dU[:, 1]

        D = p[:D]

        mul!(du, Lap, u)

        @. du = D * du + u * (1.0 - u)

        return nothing
    end,
)