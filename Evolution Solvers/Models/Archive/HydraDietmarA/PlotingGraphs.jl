using JLD2
# using GLMakie
using Plots


@load "SnapValue.jld2" DV κV UV


A1 = [1, 2, 3, 4, 5]
A2 = 6:11
A3 = 12:17

DV1 = DV[A1]
κV1 = κV[A1]
UV1 = UV[A1]

DV2 = DV[A2]
κV2 = κV[A2]
UV2 = UV[A2]

DV3 = DV[A3]
κV3 = κV[A3]
UV3 = UV[A3]



# fig = Figure()
# ax1 = Axis(fig[1, 1], xlabel = "Ω", ylabel = "U")
# lines!(ax1, Ω, UV1[1], label = "D=$(DV1[1]), κ=$(κV1[1])")


# plot(Ω, UV1, label = permutedims(["D=$(DV1[i])" for i in 1:length(UV1)]), xlabel = "Ω", ylabel = "U")

    
    
    cols5 = [
        "#ffd700",
        "#fa8775",
        "#ea5f94",
        "#cd34b5",
        "#0000ff"
    ]
        
    cols6 = [
            "#ffd700",
            "#ffb14e",
            "#fa8775",
            "#ea5f94",
            "#cd34b5",
            "#0000ff"
    ]
            
    cols7 = [
        "#ffd700",
        "#ffb14e",
        "#fa8775",
        "#ea5f94",
        "#cd34b5",
        "#9d02d7",
        "#0000ff"
    ]
                
    cols8 = [
        "#ffd700",
        "#ffb14e",
        "#fa8775",
        "#ea5f94",
        "#d646a5",
        "#cd34b5",
        "#9d02d7",
        "#0000ff"
    ]

Ω = range(0, 1, length=length(UV1[1]))  
    
function ResetPlot()
    p = plot(xlabel="x", ylabel="U(x)",     grid = true,
    minorgrid = true,
    gridalpha = 0.2,
    minorgridalpha = 0.1,
    minorticks = 4,
    guidefont = font(20),   # było ~10 → 2x większe
    tickfont = font(16),    # liczby na osiach
    legendfont = font(14)   # legenda
    )
    return p
end

function downsample(Ω, UV, m)
    N = length(Ω)
    idx = round.(Int, range(1, N, length=m))
    return Ω[idx], [u[idx] for u in UV]
end


function PlotFamily(Ω, UV; labels, colors, filename,
                    m=300, reverse_series=false,
                    lw=2, alpha=0.8)

    if m !== nothing
        Ω, UV = downsample(Ω, UV, m)
    end

    if reverse_series
        UV = reverse(UV)
        labels = reverse(labels)
    end

    p = ResetPlot()

    for i in 1:length(UV)
        plot!(p, Ω, UV[i],
            label = labels[i],
            linecolor = colors[i],
            lw = lw,
            alpha = alpha
        )
    end

    savefig(p, filename)
    display(p)
end


PlotFamily(
    Ω, UV1,
    labels = ["D=$(DV1[i])" for i in 1:length(UV1)],
    colors = cols6,
    filename = "Kappa10.pdf"
)

PlotFamily(
    Ω, UV2,
    labels = ["D=$(DV2[i])" for i in 1:length(UV2)],
    colors = cols6,
    filename = "Kappa3.pdf",
    reverse_series = true
)

PlotFamily(
    Ω, UV3,
    labels = ["κ=$(κV3[i])" for i in 1:length(UV3)],
    colors = cols7,
    filename = "Diff3.pdf",
    reverse_series = true
)
