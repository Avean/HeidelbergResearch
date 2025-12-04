using GLMakie

fig = Figure()
ax = Axis(fig[1, 1])

x = 0:0.01:2π
y = Observable(sin.(x))      # y jest Observable

lines!(ax, x, y)

fig

record(fig, "animacja.mp4", 1:200) do i
    y[] = sin.(x .+ 0.05*i)  # aktualizujemy tylko y
end