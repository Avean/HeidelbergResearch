using PyPlot

Wid = 4;
Col = :red;
# p = plot(1:10,1:10)

plot(; xlim=(0, 1), ylim=(-0.1, 1.1), grid=true, minorgrid=true)
plot!([0.2, 0.8], [0, 0], color=Col, linewidth=Wid)
plot!([0.2, 0.8], [1, 1], color=Col, linewidth=Wid)
plot!([0.2, 0.2], [0, 1], color=Col, linewidth=Wid)
plot!([0.8, 0.8], [0, 1], color=:gray, linewidth=Wid, linestyle=(0, [5, 2, 1, 2]))
plot!(legend=false)
