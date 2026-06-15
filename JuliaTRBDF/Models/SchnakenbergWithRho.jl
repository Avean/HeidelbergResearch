# models/SchnakenbergWithRho.jl

RDModel(
    id = :schnakenberg_with_rho,

    display_name = "Schnakenberg with rho(x)",

    variables = (:u, :v),

    parameters = (
        Du = 1e-4,
        Dv = 1e-2,
        a  = 0.1,
        b  = 0.9,
        ρ0 = 1.0,
        ρ1 = 0.5,
    ),

    initial = function (U, x, p)
        Random.seed!(5)

        u0 = p.a + p.b
        v0 = p.b / u0^2

        U.u .= u0 .+ 0.01 .* randn(length(x))
        U.v .= v0 .+ 0.01 .* randn(length(x))

        return nothing
    end,

    reaction = function (F, U, x, p, t)
        ρ = @. p.ρ0 + p.ρ1 * x

        @. F.u = ρ * (p.a - U.u + U.u^2 * U.v)
        @. F.v = p.b - U.u^2 * U.v

        return nothing
    end,

    diffusion = (
        u = :Du,
        v = :Dv,
    ),

    spatial_profiles = (
        ρ = (x, p) -> (@. p.ρ0 + p.ρ1 * x),
    ),
)