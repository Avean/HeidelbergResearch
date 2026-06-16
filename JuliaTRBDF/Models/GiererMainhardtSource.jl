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

τ0 = 1e8

RDModel(
    id = :gierer_meinhardt,

    display_name = "Gierer-Meinhardt system",

    variables = (:u, :v, :sd),

    parameters = (
        
        τ = τ0,
        
        Du = 1e-2,
        Dv = 1e1,
        Dsd = 1e1/τ0,

        a = 1.5,
        b = 2.0,

        μu = 0.5,
        μv = 1.0,

        pu = 0.0,
        pv = 0.0,

        ρ0 = 1.0,
        ρ1 = 0.5,
    ),

    initial = function (U, x, p)
        

        u0 = 1.0
        v0 = 2.0
        sd0 = 1.0

        U.u .= u0 .+ 0.01 .* randn(length(x))
        U.v .= v0 .+ 0.01 .* randn(length(x))
        U.sd .= sd0 .+ 0.01 .* randn(length(x))


        sd_base = @. p.ρ0 - p.ρ1 / 2 + p.ρ1 * x
        H = div(length(x), 2)

        # profile = :FootHead
        profile = :HeadFoot
        # profile = :FootHeadReverse
        # profile = :HeadFootReverse

        if profile == :FootHead
            U.sd .= sd_base

        elseif profile == :HeadFoot
            U.sd .= [sd_base[(H + 1):end]; sd_base[1:H]]

        elseif profile == :FootHeadReverse
            U.sd .= [sd_base[1:H]; reverse(sd_base[(H + 1):end])]

        elseif profile == :HeadFootReverse
            U.sd .= [reverse(sd_base[1:H]); sd_base[(H + 1):end]]
        end

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
    raw"\partial_t u = D_u \partial_{xx} u + a \cdot s \frac{u^2}{v + 1} - \mu_u u ",
    raw"\partial_t v = D_v \partial_{xx} v + b u^2 - \mu_v v",
    raw"\partial_t s = D_{s} \partial_{xx} s + \frac{1}{\tau}(u - s)",
    ),

)