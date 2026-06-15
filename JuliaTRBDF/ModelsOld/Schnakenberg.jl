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
# ============================================================

ModelSpec(
    id = :schnakenberg,

    display_name = "Schnakenberg system",

    nvars = 2,

    varnames = ["u", "v"],

    default_params = Dict(
        :Du => 1e-4,
        :Dv => 1e-2,
        :a  => 0.1,
        :b  => 0.9,
    ),

    initialize! = function (U, x, p)
        Random.seed!(2)

        N = length(x)

        a = p[:a]
        b = p[:b]

        u0 = a + b
        v0 = b / (a + b)^2

        U[:, 1] .= u0
        U[:, 2] .= v0

        U[:, 1] .+= 0.01 .* randn(N)
        U[:, 2] .+= 0.01 .* randn(N)

        return nothing
    end,

    rhs! = function (dU, U, Lap, x, p, t)
        u = @view U[:, 1]
        v = @view U[:, 2]

        du = @view dU[:, 1]
        dv = @view dU[:, 2]

        Du = p[:Du]
        Dv = p[:Dv]
        a  = p[:a]
        b  = p[:b]

        mul!(du, Lap, u)
        mul!(dv, Lap, v)

        @. du = Du * du + a - u + u^2 * v
        @. dv = Dv * dv + b - u^2 * v

        return nothing
    end,
)