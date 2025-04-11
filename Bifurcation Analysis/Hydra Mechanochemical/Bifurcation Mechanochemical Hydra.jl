using Revise 
using BifurcationKit 
using LinearAlgebra
using Plots

import LinearAlgebra: norm
const BK = BifurcationKit
using Statistics

includet("../../Evolution Solvers/Models/HydraDietmar/HydraDietmarModules.jl")

includet("../../Evolution Solvers/FillFunctions.jl")
includet("../../Evolution Solvers/DiffusionMatrices.jl")

using ..LaplaceDiscretisation
using ..FillMatrix
using ..SimParam

println("Hydra Dietmar Model - bifurcation analysis")


N = 200
L = 1
dx = L/(N-1)
X = LinRange(0, L, N) # collect(LinRange(0, L, N))

LapDis = ([-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560]);
GradDis = ([1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280]);


# we define a Bifurcation Problem
D_fixed = 1.0 / 20.0 * 1/(20*pi^2);
kappa_0 = 1+4*pi^2*D_fixed;

### Choose different for differnt approaches
# sol0 = ones(N) * kappa_0;
sol0 = [kappa_0; zeros(SimParam.SicCosNodes *2)];


par_ks = (kappa = kappa_0, Dcoef = D_fixed);

DiffMatrix = D_fixed ./ dx^2 .* FiniteDiffercePeriodic(LapDis, N);
GradMatrix = 1.0 ./ dx .* FiniteDiffercePeriodic(GradDis, N);

BasisTrun = SC.P.Trunc .* sqrt(SimParam.N-1);
BasisFull = SC.P.Full .* sqrt(SimParam.N-1);
Eigenvalues = Eig.SC;

W = ones(SimParam.N) ./ (SimParam.N - 1.0);
W[1] = W[1] / 2.0;
W[end] = W[end] / 2.0;

function F_RHS(u, par)
	(;kappa, Dcoef) = par
	return DiffMatrix*u .-u .+ kappa .*  exp.(u) ./ mean(exp.(u)) .- kappa
end

function F_SinCos(Ru, par)
	(;kappa, Dcoef) = par
	w = 
	return  - D_fixed .* Eigenvalues .* Ru .- Ru .+ kappa .*  (BasisTrun * exp.(BasisFull' * Ru) ./ (SimParam.N - 1) ) ./ (dot(exp.(BasisFull' * Ru),W)) - [kappa ; zeros(SimParam.SicCosNodes *2)];
end

function F_J_Discrete(u, par)
	(;kappa, Dcoef) = par
	return D_fixed ./ 2.0 .* mean((GradMatrix*u).^2) .+ mean(u.^2) .- kappa .*  log(mean( exp.(u))) 
end



### Choose one of the following options

# F_discr = F_RHS;
# F_discr = F_J_Discrete;
F_discr = F_SinCos;

#######


optnewton = NewtonPar(tol = 1e-11, verbose = true)
prob = BK.BifurcationProblem(F_discr, sol0, par_ks, (@optic _.kappa),
	# record_from_solution = (x, p; k...) -> (s = sum(x), u2 = x[3], nrm = norm(x)),
	record_from_solution = (x, p; k...) -> (n2 = norm(x), n8 = norm(x, 8)),
	# function to plot the solution
	plot_solution = (x, p; k...) -> plot!(x; ylabel="solution", label="", k...))
# sol = @time BK.solve( prob, Newton(), optnewton)

optcont = ContinuationPar(dsmin = 0.0001, dsmax = 0.01, ds = 0.01, p_min = 0.0, p_max = 1.6,
						  newton_options = NewtonPar(max_iterations = 30, tol = 1e-8),
						  max_steps = 300, plot_every_step = 40, n_inversion=16, nev=2*N) # , newton_options = NewtonPar(max_iterations = 10, tol = 1e-9))



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

# B =  readline()


# x = 1 .+ (1:10).^2 .* 4.0.*pi.^2.0.*D_fixed
# y =zeros(10)

# scatter!(x,y)