using Revise, BifurcationKit, LinearAlgebra, Plots
import LinearAlgebra: norm
const BK = BifurcationKit
using Statistics

includet("HydraDietmarModules.jl")
include("HydraDietmarVariables.jl")

includet("../DiffusionMatrices.jl")


using ..LaplaceDiscretisation
using ..SimParam


N = SimParam.N
L = 1
dx = L/(N-1);
X = LinRange(0, L, N) # collect(LinRange(0, L, N))


# we define a Bifurcation Problem
D_fixed = 1/(4*20*pi^2)
kappa_0 = 1+4*pi^2*D_fixed
sol0 = ones(N) * kappa_0
par_ks = (kappa = kappa_0, Dcoef = D_fixed)

DiffMatrix = D_fixed ./ dx^2 .* Lap.Per8;


# end

function F_discr(u, par)
	(;kappa, Dcoef) = par
	return DiffMatrix*u .-u .+ kappa .*  exp.(u) ./ mean(exp.(u)) .- kappa
end


##


optnewton = NewtonPar(tol = 1e-11, verbose = true)
prob = BK.BifurcationProblem(F_discr, sol0, par_ks, (@optic _.kappa),
	# record_from_solution = (x, p; k...) -> (s = sum(x), u2 = x[3], nrm = norm(x)),
	record_from_solution = (x, p; k...) -> (n2 = norm(x), n8 = norm(x, 8)),
	# function to plot the solution
	plot_solution = (x, p; k...) -> plot!(x; ylabel="solution", label="", k...))
# sol = @time BK.solve( prob, Newton(), optnewton)

optcont = ContinuationPar(dsmin = 0.0001, dsmax = 0.01, ds = 0.01, p_min = 1.0, p_max = 1.5,
						  newton_options = NewtonPar(max_iterations = 30, tol = 1e-8),
						  max_steps = 300, plot_every_step = 40, n_inversion=16, nev=N) # , newton_options = NewtonPar(max_iterations = 10, tol = 1e-9))



args = (verbosity = 0,
	plot = true,
	# callback_newton = cb, 
	halfbranch = true,
	)

function optrec(x, p, l; opt = optcont)
	level =  l
	if level <= 2
		return setproperties(opt; max_steps = 300,
			nev = 2*N, detect_loop = false)
	else
		return setproperties(opt; max_steps = 250,
			nev = 2*N, detect_loop = true)
	end
end

# automatic bifurcation diagram computation
diagram = @time bifurcationdiagram(re_make(prob, params = @set par_ks.kappa = kappa_0), PALC(),
	# very important parameter. This specifies the maximum amount of recursion
	# when computing the bifurcation diagram. It means we allow computing branches of branches
	# at most in the present case.
	2,
	optcont, bothside=true; args...
	)

p1 = plot(diagram, plotfold = true, markersize = 2, legend = :outertopleft, label = "") # legend = :outertopleft, plotstability = true, plotspecialpoints = true, putspecialptlegend = true, 
title!("#branches = $(size(diagram))")
Plots.display(p1)

##

readline()
