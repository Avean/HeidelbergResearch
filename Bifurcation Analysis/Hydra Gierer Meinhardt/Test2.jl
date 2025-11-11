using Revise
using SparseArrays
import LinearAlgebra: I, norm
using BifurcationKit
using Plots
const BK = BifurcationKit



# discretisation
N = 200
l = 6.
X = -l .+ 2l/N*(0:N-1) |> collect
h = X[2]-X[1]

# define a norm
const _weight = rand(N)
normweighted(x) = norm(_weight .* x)

# boundary condition
Δ = spdiagm(0 => -2ones(N), 1 => ones(N-1), -1 => ones(N-1) ) / h^2
L1 = -(I + Δ)^2

# functional of the problem
function R_SH!(out, u, par)
	(;λ, ν, L1) = par
	out .= L1 * u .+ λ .* u .+ ν .* u.^3 - u.^5
end

# jacobian
Jac_sp(u, par) = par.L1 + spdiagm(0 => par.λ .+ 3 .* par.ν .* u.^2 .- 5 .* u.^4)

# second derivative
d2R(u,p,dx1,dx2) = @. p.ν * 6u*dx1*dx2 - 5*4u^3*dx1*dx2

# third derivative
d3R(u,p,dx1,dx2,dx3) = @. p.ν * 6dx3*dx1*dx2 - 5*4*3u^2*dx1*dx2*dx3

# parameters associated with the equation
parSH = (λ = -0.7, ν = 2., L1 = L1)

# initial condition
sol0 = zeros(N)

# Bifurcation Problem
prob = BifurcationProblem(R_SH!, sol0, parSH, (@optic _.λ); J = Jac_sp,
	record_from_solution = (x, p; k...) -> (n2 = norm(x), nw = normweighted(x), s = sum(x), s2 = x[end ÷ 2], s4 = x[end ÷ 4], s5 = x[end ÷ 5]),
	plot_solution = (x, p;kwargs...)->(plot!(X, x; ylabel="solution", label="", kwargs...)))

opts = ContinuationPar(dsmin = 0.0001, dsmax = 0.01, ds = 0.01, p_max = 1.,
	newton_options = NewtonPar(max_iterations = 30, tol = 1e-8),
	max_steps = 300, plot_every_step = 40,
	n_inversion = 4, tol_bisection_eigenvalue = 1e-17, dsmin_bisection = 1e-7)

function cb(state; kwargs...)
	_x = get(kwargs, :z0, nothing)
	fromNewton = get(kwargs, :fromNewton, false)
	if ~fromNewton
		# if the residual is too large or if the parameter jump
		# is too big, abort continuation step
		return norm(_x.u - state.x) < 20.5 && abs(_x.p - state.p) < 0.05
	end
	true
end

args = (verbosity = 0,
	plot = true,
	callback_newton = cb, halfbranch = true,
	)

function optrec(x, p, l; opt = opts)
	level =  l
	if level <= 2
		return setproperties(opt; max_steps = 300,
			nev = N, detect_loop = false)
	else
		return setproperties(opt; max_steps = 250,
			nev = N, detect_loop = true)
	end
end    


diagram =  bifurcationdiagram(re_make(prob, params = @set parSH.λ = -0.1),
	PALC(),
	# here we specify a maximum branching level of 4
	4, optrec; args...)


plot(diagram;  plotfold = false,  
	markersize = 2, putspecialptlegend = false, xlims=(-1,1), label = "")
title!("#branches = $(size(diagram))")