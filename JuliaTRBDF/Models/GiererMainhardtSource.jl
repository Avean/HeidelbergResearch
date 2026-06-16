# models/GiererMeinhardt.jl

# ============================================================
# Classical Gierer-Meinhardt reaction-diffusion system
# ============================================================
#
# Model:
#
#     u_t = Du u_xx + a - b u + u^2 / (v+1)
#     v_t = Dv v_xx + u^2 - v
#
# Here:
#
#     u = activator
#     v = inhibitor
#
# Usually Du << Dv.
#
# ============================================================

RDModel(
    id = :gierer_meinhardt,

    display_name = "Gierer-Meinhardt system",

    variables = (:u, :v, :sd),

    parameters = (
        Du = 1e-4,
        Dv = 1e-1,
        Dsd = 1e-1,

        a = 1.5,
        b = 2.0,
        τ = 1.0,

        μu = 0.5,
        μv = 1.0,

        pu = 0.0,
        pv = 0.0,

        ρ0 = 1.0,
        ρ1 = 0.5,
    ),

    initial = function (U, x, p)
        Random.seed!(6)

        u0 = (p.a + 1.0) / p.b
        v0 = u0^2
        sd0 = u0

        U.u .= u0 .+ 0.01 .* randn(length(x))
        U.v .= v0 .+ 0.01 .* randn(length(x))
        U.sd .= sd0 .+ 0.01 .* randn(length(x))

        return nothing
    end,

    reaction = function (F, U, x, p, t)
        @. F.u = p.a * U.u^2 / (U.v + 1.0) - p.μu * U.u + p.pu
        @. F.v = p.b * U.u^2 - p.μv * U.v + p.pv
        @. F.sd = (U.u - U.sd) / p.τ

        return nothing
    end,

    diffusion = (
        u = :Du,
        v = :Dv,
        sd = :Dsd,
    ),

    latex_equations = (
    raw"\partial_t u = D_u \partial_{xx} u + a \frac{u^2}{v + 1} - \mu_u u + p_u",
    raw"\partial_t v = D_v \partial_{xx} v + b u^2 - \mu_v v + p_v",
    ),

)