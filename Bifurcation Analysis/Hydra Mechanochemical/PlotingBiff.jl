D = diagram.child[1]

x = D.γ.branch.param
y = D.γ.branch.n2

p = plot(x,y, legend = false, color = "lightgray")

x0 = x[1];
x1 = x[48:160]
y1 = y[48:160]


for D in diagram.child
    display(1)
    x = D.γ.branch.param
    y = D.γ.branch.n2
    plot!(x,y, color = "lightgray", linewidth = 2.0)
end

xlims!(0.3,1.6)


plot!(x1,y1, color = "black", linewidth = 3.0)
plot!([0, x0],[0.0, 0.0], color = "black", linewidth = 3.0)
plot!([x0, 1.6],[0.0, 0.0], color = "lightgray", linewidth = 2.0)



display(p)