using GLMakie
using Makie

fig = Figure()

ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])

title1 = Observable("Sinus")
title2 = Observable("Punkty")

h1 = lines!(ax1, 1:10, sin.(1:10), label = lift(t -> t, title1))
h2 = scatter!(ax2, 1:10, rand(10), label = lift(t -> t, title2))

AxisLegend(ax1)
AxisLegend(ax2)

fig