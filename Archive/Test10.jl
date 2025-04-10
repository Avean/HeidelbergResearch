using Revise, BifurcationKit, LinearAlgebra, Plots
const BK = BifurcationKit

N(x; a = 0.5, b = 0.01) = 1 + (x + a*x^2)/(1 + b*x^2)

function F_chan(x, p)
	(;α, β) = p
	f = similar(x)
	n = length(x)
	f[1] = x[1] - β
	f[n] = x[n] - β
	for i=2:n-1
		f[i] = (x[i-1] - 2 * x[i] + x[i+1]) * (n-1)^2 + α * N(x[i], b = β)
	end
	return f
end
n = 101
sol0 = [(i-1)*(n-i)/n^2+0.1 for i=1:n]

# set of parameters
par = (α = 3.3, β = 0.01)

optnewton = NewtonPar(tol = 1e-11, verbose = true)

prob = BifurcationProblem(F_chan, sol0, par, (@optic _.α),
	# function to plot the solution
	plot_solution = (x, p; k...) -> plot!(x; ylabel="solution", label="", k...))
sol = @time BK.solve( prob, Newton(), optnewton)


optcont = ContinuationPar(dsmin = 0.01, dsmax = 0.2, ds= 0.1, p_min = 0., p_max = 4.2,
	newton_options = NewtonPar(max_iterations = 10, tol = 1e-9))


    br = continuation(prob, PALC(), optcont; plot = true)
nothing #hide