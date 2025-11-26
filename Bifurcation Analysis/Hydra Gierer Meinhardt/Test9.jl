using GLMakie

@time F = Figure()
@time A = Axis(F[1, 1])
@time B = Slider(F[1, 2], range = 0.0:0.1:10.0, startvalue = 1.0)
@time C = Textbox(F[1,3], placeholder = "Test")
N = 1e4
@time scatter!(A, rand(Int(N)), zeros(Int(N)))

scatter!(A,Float64[],Float64[])