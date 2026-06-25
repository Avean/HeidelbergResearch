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

τ0 = 1e0

RDModel(
    id = :gierer_meinhardt,

    display_name = "Gierer-Meinhardt system",

    variables = (:u, :v, :h, :sd),

    parameters = (
        
        τ = τ0,
        
        Du = 1e-2,
        Dv = 1e0,
        Dh = 0.0,
        Dsd = 1e0/τ0,

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
        sd0 = 1.0
        h0 = 0.0

        U.u .= u0 .+ 0.01 .* randn(length(x))
        U.v .= v0 .+ 0.01 .* randn(length(x))
        U.h .= h0 
        


        sd_base = @. p.ρ0 - p.ρ1 / 2 + p.ρ1 * x
        H = div(length(x), 2)

        # profile = :Constant
        # profile = :FootHead
        # profile = :HeadFoot
        # profile = :FootHeadReverse
        # profile = :HeadFootReverse
        profile = :ZigZag


        if profile == :FootHead
            U.sd .= sd_base

        elseif profile == :HeadFoot
            U.sd .= [sd_base[(H + 1):end]; sd_base[1:H]]

        elseif profile == :FootHeadReverse
            U.sd .= [sd_base[1:H]; reverse(sd_base[(H + 1):end])]

        elseif profile == :HeadFootReverse
            U.sd .= [reverse(sd_base[1:H]); sd_base[(H + 1):end]]
        
        elseif profile == :ZigZag
            x1 = 0.2
            x2 = 0.9

            s = p.ρ1

            U.sd .= @. ifelse(
                x < x1,
                p.ρ0 + s * x,
                ifelse(
                    x < x2,
                    p.ρ0 + s * (x - x1),
                    p.ρ0 + s * (x - x2),
                )
            )
        elseif profile == :Constant
            U.sd .= sd0 .+ 0.01 .* randn(length(x))
        end

        

        return nothing
    end,

    

    reaction = function (F, U, x, p, t)
        @. F.u = p.a * U.sd * U.u^2 / (U.v + 1.0) - p.μu * U.u + p.pu
        @. F.v = p.b * U.sd * U.u^2 - p.μv * U.v + p.pv
        # @. F.h = 0.01 .* (U.u - 14.0 * U.h *(U.h^2.0 .* 4.0 - 5.2 * U.h + 2.0))
        @. F.h = 0.01 .*(U.u - U.h)
        @. F.sd = (1.0 + 5*U.h  - U.sd) / p.τ

        return nothing
    end,

    diffusion = (
        u = :Du,
        v = :Dv,
        h = :Dh,
        sd = :Dsd,
    ),

    latex_equations = (
    raw"\partial_t u = D_u \partial_{xx} u + a \cdot s \frac{u^2}{v + 1} - \mu_u u ",
    raw"\partial_t v = D_v \partial_{xx} v + b \cdot s u^2 - \mu_v v",
    raw"\gamma\partial_t h = D_h \partial_{xx} h + u - h",
    raw"\tau\partial_t s = D_{s} \partial_{xx} s + 1 + c \cdot h - s",
    ),

)