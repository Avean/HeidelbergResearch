using Plots
using Polynomials

function phase_portrait(f, g;
    u0 = -2.0,
    u1 =  2.0,
    v0 = -2.0,
    v1 =  2.0,
    Nu = 25,
    Nv = 25,
    arrow_length = 0.15,
    h = nothing,
    equilibria = Tuple{Float64,Float64,String}[]
)
    # Siatka punktów
    us = range(u0, u1, length = Nu)
    vs = range(v0, v1, length = Nv)

    # Punkty startowe strzałek
    U = Float64[]
    V = Float64[]

    # Składowe strzałek
    dU = Float64[]
    dV = Float64[]

    for u in us
        for v in vs
            fu = f(u, v)
            gv = g(u, v)

            norm = sqrt(fu^2 + gv^2)
            L = 1
            if norm > 1e-12
                push!(U, u)
                push!(V, v)
                push!(dU, arrow_length * fu / norm/L)
                push!(dV, arrow_length * gv / norm/L)
            end
        end
    end

    # Bazowy wykres: portret fazowy
    p = quiver(
        U, V;
        quiver = (dU ./ 2, dV ./ 2),
        xlabel = "u",
        ylabel = "v",
        xlims = (u0, u1),
        ylims = (v0, v1),
        linewidth = 0.5,
        aspect_ratio = :equal,
        legend = :bottomright,
        label = ""
    )

    # Dorysowanie krzywej v = h(u)
    if h !== nothing
        u_curve = range(u0, u1, length = 1000)
        v_curve = [h(u) for u in u_curve]

        plot!(
            p,
            u_curve,
            v_curve;
            color = :black,
            linewidth = 1,
            label = ""
        )
    end

    # Dorysowanie punktów stacjonarnych według typu
    if !isempty(equilibria)

        classes = [
            ("Stable",   :green),
            ("Unstable", :red),
            ("DDI",      :blue)
        ]
        # Kolory punktów stacjonarnych według statusu
        status_colors = Dict(
            "stable"   => :green,
            "unstable" => :red,
            "DDI"      => :blue
        )

        # Pamiętamy, które statusy już pojawiły się w legendzie
        labels_used = Set{String}()

        for eq in equilibria
            u_eq = eq[1]
            v_eq = eq[2]

            # Jeśli punkt ma status, bierzemy go z trzeciej współrzędnej.
            # Jeśli nie ma statusu, traktujemy go jako zwykły punkt.
            status = length(eq) >= 3 ? eq[3] : "equilibrium"

            # Kolor zależy od statusu.
            # Jeśli status nie jest znany, używamy czarnego.
            point_color = get(status_colors, status, :black)

            # Etykieta w legendzie ma pojawić się tylko raz dla danego statusu.
            point_label =
                if status in labels_used
                    ""
                else
                    push!(labels_used, status)
                    status
                end

            scatter!(
                p,
                [u_eq],
                [v_eq];
                color = point_color,
                markersize = 5,
                markerstrokecolor = :black,
                label = point_label
            )
        end
    end

    return p
end


μu = 0.5
μv = 1.0
a  = 1.5
b  = 2.0
pu = 0.0
pv = 0.0
K = 1.0
cu = 0.01

f(u, v) = -μu * u + a * (u^2+ cu) / (K + v) + pu
g(u, v) = -μv * v + b * u^2

# Przykładowa krzywa: nullklina g(u,v)=0, czyli v = (b*u^2 + pv)/μv
h(u) = (b * u.^2 .+ pv) ./ μv

z(u) = a*(u^2 + cu) -  μu * u * (K + h(u))
p = fit([0.0, 1.0, 2.0, 3.0], z.([0.0, 1.0, 2.0, 3.0]), 3)
PP = real(roots(p))

println(roots(p))


# Punkty stacjonarne z etykietami:
# (u, v, "stable")
# (u, v, "unstable")
# (u, v, "DDI")
# equilibria = [
#     (0.0, 0.0, "Stable"),
#     (0.5, 0.5, "Unstable"),
#     (1.0, 2.0, "DDI")
# ]
equilibria = [
    (PP[1], h(PP[1]), "stable"),
    (PP[2], h(PP[2]), "unstable"),
    (PP[3], h(PP[3]), "DDI")
]


p = phase_portrait(f, g;
    u0 = -0.1,
    u1 = 2.0,
    v0 = -0.1,
    v1 = 2.7,
    Nu = 20,
    Nv = 20,
    h = h,
    equilibria = equilibria
)

hf(u) = (a / μu) * (u^2 + cu) / u - K
u_curve = range(1e-4, 2.0, length = 1000)

plot!(
    p,
    u_curve,
    hf.(u_curve);
    color = :blue,
    linewidth = 1,
    label = ""
)

savefig(p, joinpath(@__DIR__, "phase_portrait3.pdf"))