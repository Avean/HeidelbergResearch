# models/ThreeVariableTest.jl

# ============================================================
# Three-variable test reaction-diffusion system
# ============================================================
#
# A simple test model with three variables:
#
#     u_t = Du u_xx + a - u + u^2 v + gamma * (s - u)
#     v_t = Dv v_xx + b - u^2 v
#     s_t = Ds s_xx + eps * (u - s)
#
# The third variable s acts as a slowly relaxing source/background field.
#
# ============================================================

RDModel(
    id = :three_variable_test,

    display_name = "Three-variable test system",

    variables = (:u, :v, :s),

    parameters = (
        Du    = 1e-4,
        Dv    = 1e-2,
        Ds    = 5e-3,
        a     = 0.1,
        b     = 0.9,
        gamma = 0.2,
        eps   = 0.05,
    ),

    initial = function (U, x, p)
        Random.seed!(3)

        u0 = p.a + p.b
        v0 = p.b / u0^2
        s0 = u0

        U.u .= u0 .+ 0.01 .* randn(length(x))
        U.v .= v0 .+ 0.01 .* randn(length(x))
        U.s .= s0 .+ 0.01 .* cos.(2π .* x)

        return nothing
    end,

    reaction = function (F, U, x, p, t)
        @. F.u = p.a - U.u + U.u^2 * U.v + p.gamma * (U.s - U.u)
        @. F.v = p.b - U.u^2 * U.v
        @. F.s = p.eps * (U.u - U.s)

        return nothing
    end,

    diffusion = (
        u = :Du,
        v = :Dv,
        s = :Ds,
    ),
)