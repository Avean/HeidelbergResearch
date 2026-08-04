# models/StemCells/StemCells.jl

# ============================================================
# Stem-cell reaction-diffusion system
# ============================================================
#
# Model:
#
#   Q_t = dQ Q_xx
#         - r0 Q² / (K + Q²)
#         + 2 b0 p0 A(Q + A) / [(1 + beta Q)(H + Q + A)]
#
#   A_t = dA A_xx
#         + gamma r0 Q² / (K + Q²)
#         - p0 A(Q + A) / (H + Q + A)
#
# Here:
#
#   Q = quiescent stem-cell population
#   A = active stem-cell population
#
# ============================================================

RDModel(
    id = :stem_cells,

    display_name = "Stem-cell system",

    variables = (:Q, :A),

    parameters = (
        DQ = 0.001,
        DA = 0.9,

        r0 = 9.999943,
        K = 489.290158,

        p0 = 999.998129,
        H = 999.998179,

        beta = 0.000669497,
        b0 = 0.715846,
        gamma = 0.723547,
    ),

    initial = function (U, x, p)
        # Approximate positive homogeneous equilibrium
        Q0 = 55.0
        A0 = 0.11

        U.Q .= Q0 .* (1.0 .+ 0.001 .* randn(length(x)))
        U.A .= A0 .* (1.0 .+ 0.001 .* randn(length(x)))

        return nothing
    end,

    reaction = function (F, U, x, p, t)
        @. F.Q =
            -p.r0 * U.Q^2 / (p.K + U.Q^2) +
            2.0 * p.b0 * p.p0 * (U.Q + U.A) * U.A /
            ((1.0 + p.beta * U.Q) * (p.H + U.Q + U.A))

        @. F.A =
            p.gamma * p.r0 * U.Q^2 / (p.K + U.Q^2) -
            p.p0 * (U.Q + U.A) * U.A /
            (p.H + U.Q + U.A)

        return nothing
    end,

    diffusion = (
        Q = :DQ,
        A = :DA,
    ),

    latex_equations = (
        raw"\partial_t Q = D_Q \partial_{xx} Q
        - r_0 \frac{Q^2}{K+Q^2}
        + \frac{2b_0p_0(Q+A)A}{(1+\beta Q)(H+Q+A)}",

        raw"\partial_t A = D_A \partial_{xx} A
        + \gamma r_0 \frac{Q^2}{K+Q^2}
        - p_0\frac{(Q+A)A}{H+Q+A}",
    ),
)
