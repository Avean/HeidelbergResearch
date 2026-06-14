# models/SchnakenbergSource.jl

# ============================================================
# Schnakenberg system with source density
# ============================================================
#
# Model:
#
#     u_t   = Du u_xx   + rho * (a - u + u^2 v)
#     v_t   = Dv v_xx   + b - u^2 v
#     rho_t = Drho rho_xx + eps * (rho0 - rho)
#
# Variables:
#
#     u   - activator
#     v   - inhibitor
#     rho - source density
#
# This is a simple three-equation demo model. The source density rho
# modulates the production/reaction term of the activator.
#
# ============================================================

ModelSpec(
    id = :schnakenberg_source,

    display_name = "Schnakenberg + source density",

    nvars = 3,

    varnames = ["u", "v", "rho"],

    default_params = Dict(
        :Du      => 1e-4,
        :Dv      => 1e-2,
        :Drho    => 5e-4,
        :a       => 0.1,
        :b       => 0.9,
        :rho0    => 1.0,
        :rho_amp => 0.3,
        :eps     => 0.05,
    ),

    initialize! = function (U, x, p)
        Random.seed!(3)

        N = length(x)

        a       = p[:a]
        b       = p[:b]
        rho0    = p[:rho0]
        rho_amp = p[:rho_amp]

        u0 = a + b
        v0 = b / (a + b)^2

        U[:, 1] .= u0
        U[:, 2] .= v0

        # Smooth nonconstant source density profile.
        # The derivative of cos(pi*x) vanishes at x = 0 and x = 1,
        # so it is compatible with homogeneous Neumann boundary conditions.
        U[:, 3] .= @. rho0 + rho_amp * cos(pi * x)

        U[:, 1] .+= 0.01 .* randn(N)
        U[:, 2] .+= 0.01 .* randn(N)
        U[:, 3] .+= 0.005 .* randn(N)

        return nothing
    end,

    rhs! = function (dU, U, Lap, x, p, t)
        u   = @view U[:, 1]
        v   = @view U[:, 2]
        rho = @view U[:, 3]

        du   = @view dU[:, 1]
        dv   = @view dU[:, 2]
        drho = @view dU[:, 3]

        Du   = p[:Du]
        Dv   = p[:Dv]
        Drho = p[:Drho]

        a    = p[:a]
        b    = p[:b]
        rho0 = p[:rho0]
        eps  = p[:eps]

        mul!(du,   Lap, u)
        mul!(dv,   Lap, v)
        mul!(drho, Lap, rho)

        @. du   = Du   * du   + rho * (a - u + u^2 * v)
        @. dv   = Dv   * dv   + b - u^2 * v
        @. drho = Drho * drho + eps * (rho0 - rho)

        return nothing
    end,
)