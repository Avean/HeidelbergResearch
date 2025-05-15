D = diagram.child[1]

x = D.γ.branch.param
y = D.γ.branch.n2

plot(x,y)





x0 = 1 + 4*pi^2*D_fixed;
xbiff = x[argmin(x)]
ybiff = y[argmin(x)]

x1 = x[argmin(x):end]
y1 = y[argmin(x):end]

# x1 = x[argmax(x):argmin(x)]
# y1 = y[argmax(x):argmin(x)]


Red = RGBA(0.917, 0.353, 0.18, 0.3);
Red = RGB(0.98, 0.72, 0.70);
Blue = "#0072BD"

p = plot(x,y, legend = false, color = Red)

for D in diagram.child
    x2 = D.γ.branch.param
    y2 = D.γ.branch.n2

    xbiff = [xbiff; x2[argmin(x2)]]
    ybiff = [ybiff; y2[argmin(x2)]]    
    plot!(x2,y2, color = Red, linewidth = 2.0)
end




plot!(x1,y1, color = Blue, linewidth = 3.0)
plot!([x0, 12.6],[0.0, 0.0], color = Red, linewidth = 2.0)
plot!([0, x0],[0.0, 0.0], color = Blue, linewidth = 3.0)


#### Biffuracion points

x = 1 .+  4.0.*pi^2.0.*(1:10).^2.0.*D_fixed;
y = zeros(length(x));

x = [x; xbiff]
y = [y; ybiff]

kolor_szary = RGB(0.7, 0.7, 0.7)         # Jasny szary (środek)
kolor_czarny = :black                    # Czarna obwódka

# Wykres scatter
scatter!(
    x, y,
    markercolor=kolor_szary,             # Wypełnienie markerów
    markerstrokecolor=kolor_czarny,      # Obramowanie markerów
    markersize=3.5,
    markerstrokewidth=1,
    legend=false,
)

#Diffusion coefficent 0.0004 and 15 Fourier nodes
# plot!(p, tickfont = 16, grid = true, minorgrid = true, 
#         xticks = (0.5:0.5:2.5, [0.5, "", 1.5, "", 2.5,]),
#         yticks = (0:5, [""]),        
#         xlims = (0.3, 2.6),
#         ylims = (-0.5, 5.0),
#         )


# Diffusion coefficent 0.02, and 3 Fourier nodes
plot!(p, tickfont = 16, grid = true, minorgrid = true, 
        xticks = (0:2.5:10.0, [0.0, "", 5.0, "", 10.0]),
        yticks = (0:7, [""]),        
        xlims = (0.0, 10.5),
        ylims = (-0.5, 6.0),
        )

display(p)