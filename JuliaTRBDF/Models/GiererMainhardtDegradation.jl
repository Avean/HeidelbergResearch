# models/GiererMeinhardt.jl

# ============================================================
# Classical Gierer-Meinhardt reaction-diffusion system with degradation
# ============================================================
#
# Model:
#
#     u_t = Du u_xx + a - b u + u^2 - u(v+1) + ρ_a
#     v_t = Dv v_xx + u^2 - v + ρ_v
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

    variables = (:u, :v),

    parameters = (
        Du = 1e-2,
        Dv = 1e2,

        a = 0.75,
        b = 0.5,

        μu = 0.5,
        μv = 2.0,

        pu = 0.0,
        pv = 0.0,

        ρ0 = 1.0,
        ρ1 = 0.5,
    ),

    initial = function (U, x, p)
        # Random.seed!(6)

        u0 = 2.0
        v0 = 1.0

        U.u .= u0 .+ 0.01 .* randn(length(x))
        U.v .= v0 .+ 0.01 .* randn(length(x))

        return nothing
    end,

    reaction = function (F, U, x, p, t)
        @. F.u = p.a * U.u^2 - U.u * U.v - p.μu * U.u + p.pu 
        @. F.v = p.b * U.u^2 - p.μv * U.v + p.pv

        return nothing
    end,

    diffusion = (
        u = :Du,
        v = :Dv,
    ),

    latex_equations = (
    raw"\partial_t u = D_u \partial_{xx} u + a u^2 - uv - \mu_u u ",
    raw"\partial_t v = D_v \partial_{xx} v + b u^2 - \mu_v v ",
    ),

)