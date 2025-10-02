X1 = getfield(W1, Fields[1]);
X2 = getfield(W2, Fields[1]);
X3 = getfield(W3, Fields[1]);

p= plot(layout = (1,3),
    size = (1200, 300),
    bottom_margin=5mm,
)

plot!(p,subplot = 1,SimParam.x,X1, 
        color = :black,
        linewidth = 2.0,
        # size = (100,100),
        )

plot!(p, subplot = 1, legend = false,
    grid = true, minorgrid = true,
    # title = "D = $D, κ = $κ",
    xlims = (-0.05, 1.05),
    xticks = (0.0:0.2:1.0, [0.0, "","", "", "", 1.0]),
    yticks = (0:5:20, ["0.0", "", "", "", "20.0"]),
    tickfont = 16,
    ylims = (0,20),
    )

plot!(p, subplot = 2,SimParam.x,X2, 
        color = :black,
        linewidth = 2.0,)

plot!(p, subplot = 2, legend = false,
    grid = true, minorgrid = true,
    # title = "D = $D, κ = $κ",
    xlims = (-0.05, 1.05),
    xticks = (0.0:0.2:1.0, [0.0, "","", "", "", 1.0]),
    yticks = (0:10:50, ["0.0", "", "", "", "","50.0"]),
    tickfont = 16,
    ylims = (0,50),
    )

plot!(p, subplot = 3,SimParam.x,X3, 
    color = :black,
    linewidth = 2.0,)

plot!(p, subplot = 3,legend = false,
grid = true, minorgrid = true,
# title = "D = $D, κ = $κ",
xlims = (-0.05, 1.05),
xticks = (0.0:0.2:1.0, [0.0, "","", "", "", 1.0]),
yticks = (0:50:150, ["0.0", "", "", "150.0"]),
tickfont = 16,
ylims = (0,150),
)

# 
# tight_layout()