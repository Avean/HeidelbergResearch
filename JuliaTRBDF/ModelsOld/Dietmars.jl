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
    id = :Dietmar,

    display_name = "Dietmar Mechanochemical",

    nvars = 1,

    varnames = ["u"],

    default_params = Dict(
        :D => 0.01,
        :κ => 2,
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
        κ = p[:κ]

        mul!(du, Lap, u)

        # Approximation of ∫_0^1 u(x) dx by the trapezoidal rule.
        IEu = (x[2]-x[1]) * (sum(exp.(u)) - 0.5 * (exp(u[1]) + exp(u[end])))

        # Safety guard against division by zero.
        IEu = max(IEu, eps(Float64))

        @. du = D * du + κ * exp(u)/IEu - u

        return nothing
    end,
)