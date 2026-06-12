using Plots

function phase_portrait(f, g;
    u0 = -2.0,
    u1 =  2.0,
    v0 = -2.0,
    v1 =  2.0,
    Nu = 25,
    Nv = 25,
    arrow_length = 0.15,
    h = nothing,
    equilibria = Tuple{Float64,Float64,String}[],
    title = "Phase portrait"
)
    us = range(u0, u1, length = Nu)
    vs = range(v0, v1, length = Nv)

    U  = Float64[]
    V  = Float64[]
    dU = Float64[]
    dV = Float64[]

    for u in us
        for v in vs
            fu = f(u, v)
            gv = g(u, v)

            n = sqrt(fu^2 + gv^2)

            if n > 1e-12
                push!(U, u)
                push!(V, v)
                push!(dU, arrow_length * fu / n)
                push!(dV, arrow_length * gv / n)
            end
        end
    end

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
        title = title,
        label = ""
    )

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

    if !isempty(equilibria)
        classes = [
            ("Stable",   :green),
            ("Unstable", :red),
            ("DDI",      :blue)
        ]

        for (class_name, class_color) in classes
            pts = [(u, v) for (u, v, status) in equilibria if status == class_name]

            if !isempty(pts)
                u_eq = [pt[1] for pt in pts]
                v_eq = [pt[2] for pt in pts]

                scatter!(
                    p,
                    u_eq,
                    v_eq;
                    color = class_color,
                    markersize = 4,
                    markerstrokecolor = :black,
                    label = class_name
                )
            end
        end
    end

    return p
end