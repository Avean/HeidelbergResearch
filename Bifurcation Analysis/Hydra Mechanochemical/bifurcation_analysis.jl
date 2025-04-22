using Revise, ForwardDiff
using BifurcationKit, LinearAlgebra, Plots, Statistics, SparseArrays
const BK = BifurcationKit

###############################################################################################################################

# Define the discretisation of our Problem
function Laplacian1D(Nx, lx)
	hx = 2 * lx / Nx
	D2x = spdiagm(0 => -2 * ones(Nx), 1 => ones(Nx-1), -1 => ones(Nx-1)) / hx^2
	D2x[end, 1] = 1 / hx^2
	D2x[1, end] = 1 / hx^2

	D2xsp = sparse(D2x)
	return D2xsp
end

ϕ(u, κ)  = - u + κ*exp.(u)/mean(exp.(u))
delta(n) = n == 0
dϕ(u, κ, i, k) = (-1 + κ*exp(u[i]) / mean(exp.(u))) * delta(i-k) - κ*exp(u[i]) .* exp(u[k]) / (mean(exp.(u))*sum(exp.(u)))

function Fmit!(f, u, p)
	mul!(f, p.Δ, u)
	f .= p.diffcoef*f .+ ϕ(u, p.κ)
	return f
end

function JFmit(x,p)
	J = p.Δ
	dNonlin = zeros(length(x), length(x))
	for i in 1:length(x)
		for k in 1:length(x)
			dNonlin[i, k] = dϕ(x, p.κ, i, k)
		end
	end
	dNonlin = sparse(dNonlin)
	return p.diffcoef*J + dNonlin
end

function F_discr!(f, u, p)
	mul!(f, p.Δ, u)
	f .= p.diffcoef*f .- u .+ p.κ * exp.(u) ./ mean(exp.(u))# - p.κ*ones(length(u))
	return f
end

###########################################################################################################################
# Bifurcation Problem
# Parameters:
# DiffCoef = 1/(40*pi^2);
DiffCoef = 1/ (20*pi^2);
# DiffCoef = 1/ (2*pi^2);
kappa_0 = 1 + 4 * pi^2 * DiffCoef;
Nx = 100; lx = 0.5;
Δ = Laplacian1D(Nx, lx);
par_mit = (diffcoef = DiffCoef, κ = kappa_0, Δ = Δ);
# Define the L2 norm with weight.
# Choose the weighted norm in order to break symmetries and see all branches.
w = (lx .+ LinRange(-lx, lx, Nx)) |> vec
w_one = ones(Nx) |> vec
norm2(x) = norm(x .* w_one) / sqrt(length(x))
norm2_weighted(x) = norm(x .* w) / sqrt(length(x))
# initial guess f for newton
sol0 = zeros(Nx) |> vec

nev_N = Nx - 2
int_kappa = [0., 1 + 9*4*pi^2*DiffCoef + kappa_0/4.]

prob = BifurcationProblem(Fmit!, sol0, par_mit, (@optic _.κ),; J=JFmit,
  record_from_solution = (x, p; k...) -> (nrm = norm2(x), nw = norm2_weighted(x), n∞ = norminf(x), sol=x),
  plot_solution = (x, p; k...) -> plot!(x ; k...))


# eigensolver
eigls = EigArpack()

# options for Newton solver, we pass the eigen solver
opt_newton = BK.NewtonPar(tol = 1e-8, eigsolver = eigls, max_iterations = 20)

# options for continuation
opts_br = ContinuationPar(p_max = int_kappa[2], p_min = int_kappa[1],
	# for a good looking curve
	dsmin = 0.001, dsmax = 0.05, ds = 0.01,
	# number of eigenvalues to compute
	nev = nev_N,
	newton_options = (@set opt_newton.verbose = false),
	tol_stability = 1e-6,
	# detect codim 1 bifurcations
	detect_bifurcation = 3,
	# Optional: bisection options for locating bifurcations
	n_inversion = 4, dsmin_bisection = 1e-7, max_bisection_steps = 25)

# optional arguments for continuation
kwargsC = (verbosity = 0, plot = true, normC = norminf)

br = continuation(prob, PALC(), opts_br, bothside=true; kwargsC...)

all_branches = Array{Vector{Branch}}(undef, length(br.specialpoint)-2)
for (index, point) in pairs(br.specialpoint)
	@show index
	@show point.param
	if index == 1 || index == length(br.specialpoint)
		# Those are only endpoints of interval for κ and not bifurcations points, so skip them
		continue
	end
	branches = continuation(br, index,
		setproperties(opts_br; ds = 0.01), bothside=true ;
		alg = PALC(),
		kwargsC...,
		nev = nev_N,
	)
	if branches == Branch[]
		@show "No branches found for index $(index), try smaller ds"
		branches = continuation(br, index,
		setproperties(opts_br; ds = 0.001),
		bothside=true ;
		alg = PALC(),
		kwargsC...,
		nev = nev_N,
		)
	end
	all_branches[index-1] = branches
end

p1 = plot(br)
for branch in all_branches
	p1 = plot!(branch...)
end
display(p1)

# Use weighted norm for plotting to see different branches due to the symmetry breaking
p2 = plot(br; vars = (:param, :nw), title = "Weighted norm")
for branch in all_branches
	p2 = plot!(branch...; vars = (:param, :nw))
end
plot(p1, p2)

######################################################################################################################
# To have a clean bifurcation diagram, we only plot one branch of the symmetric branches.
# Choose branch[1] for the first branch, branch[2] for the second branch...
# Maybe have to play around with the index of the branch to get the one you want. Can be different for different parameters.
# Can also be differnt for different bifurcation points.
p3 = plot(br)
for branch in all_branches
	p3 = plot!(branch[1])
end
display(p3)
# ylims!(p3, 18, 20) # Specify the y-axis limits for the plot
# xlims!(p3, 15, 20) # Specify the x-axis limits for the plot

######################################################################################################################
# One can also use bifurcationdiagram() to compute the bifurcation diagram automatically and then by hand the missing branches.
# automatic bifurcation diagram computation
diagram = @time bifurcationdiagram(prob, PALC(),
	# very important parameter. This specifies the maximum amount of recursion
	# when computing the bifurcation diagram. It means we allow computing branches of branches
	# at most in the present case.
	2,
	opts_br, bothside=true;
	kwargsC...
	)
p4 = plot(diagram, plotfold = true, markersize = 2, label = "", putspecialptlegend = false)
# ylims!(p4, 0, 5)  # Specify the y-axis limits for the plot

# If some branch in automatic bifurcation diagram is not computed, we can compute it by hand
# (1 = endpoint; 2,... = bifurcation points; length(diagram.γ.specialpoint) = endpoint):
br_missing = continuation(diagram.γ, 3,
		setproperties(opts_br; detect_bifurcation = 3, ds = 0.001, p_min = int_kappa[1], p_max = int_kappa[2]), bothside=true ;
		alg = PALC(),
		kwargsC...,
		nev = nev_N,
		)
plot!(p4, br_missing...)

# Plot the bifurcation diagram with weighted norm to see different branches	due to the symmetry breaking
p5 = plot(diagram; putspecialptlegend=false, markersize=2, vars = (:param, :nw), label = "",
		  title = "Plot with with weighted norm, #branches = $(size(diagram))")

######################################################################################################################
