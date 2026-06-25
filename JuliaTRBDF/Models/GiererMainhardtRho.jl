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

    variables = (:u, :v),

    parameters = (
        Du = 1e-2,
        Dv = 1e0,

        a = 1.5,
        b = 2.0,

        μu = 0.5,
        μv = 1.0,

        pu = 0.0,
        pv = 0.0,

        ρ0 = 1.0,
        ρ1 = 1.5,
        # ρ1 = 0.5,
    ),

    initial = function (U, x, p)

        u0 = 1.0
        v0 = 2.0

        U.u .= u0 .+ 0.01 .* randn(length(x))
        U.v .= v0 .+ 0.01 .* randn(length(x))

        return nothing
    end,

    


    reaction = function (F, U, x, p, t)
        @. F.u = p.a * p.ρ * U.u^2 / (U.v + 1.0) - p.μu * U.u + p.pu
        @. F.v = p.b * p.ρ * U.u^2 - p.μv * U.v + p.pv

        return nothing
    end,

    diffusion = (
        u = :Du,
        v = :Dv,
    ),




    spatial_profiles = (
        FootHead = (
            ρ = (x, p) -> begin
                ρx = @. p.ρ0 - p.ρ1 / 2 + p.ρ1 * x
                return ρx
            end,
        ),

        HeadFoot = (
            ρ = (x, p) -> begin
                H = div(length(x), 2)
                ρx = @. p.ρ0 - p.ρ1 / 2 + p.ρ1 * x

                return [ρx[(H + 1):end]; ρx[1:H]]
            end,
        ),

        FootHeadReverse = (
            ρ = (x, p) -> begin
                H = div(length(x), 2)
                ρx = @. p.ρ0 - p.ρ1 / 2 + p.ρ1 * x

                return [ρx[1:H]; reverse(ρx[(H + 1):end])]
            end,
        ),

        HeadFootReverse = (
            ρ = (x, p) -> begin
                H = div(length(x), 2)
                ρx = @. p.ρ0 - p.ρ1 / 2 + p.ρ1 * x

                return [reverse(ρx[(1):H]); (ρx[(H+1):end])]
            end,
        ),

        ZigZag = (
            ρ = (x, p) -> begin
                x1 = 0.2
                x2 = 0.9

                ρx = @. ifelse(
                    x < x1,
                    p.ρ0 + p.ρ1 * x,
                    ifelse(
                        x < x2,
                        p.ρ0 + p.ρ1 * (x - x1),
                        p.ρ0 + p.ρ1 * (x - x2),
                    )
                )

                return ρx
            end,
        ),
    ),

    latex_equations = (
    raw"\partial_t u = D_u \partial_{xx} u + a\cdot\rho(x) \frac{u^2}{v + 1} - \mu_u u",
    raw"\partial_t v = D_v \partial_{xx} v + b\cdot \rho(x) u^2 - \mu_v v",
    ),
)