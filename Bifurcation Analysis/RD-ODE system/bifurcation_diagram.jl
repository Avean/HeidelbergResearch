using Revise, ForwardDiff
using BifurcationKit, LinearAlgebra, Plots, Statistics, SparseArrays
const BK = BifurcationKit

###############################################################################################################################

# Define the discretisation of our Problem
# Laplacian operator in 1D with Neumann boundary conditions

f1(v, w) = v * v * w

function Fbru!(f, x, p, t = 0)
	(;pp, qq, D1, D2, l) = p
	n = div(length(x), 3)
	h2 = 1.0 / n^2
	c1 = D1 / l^2 / h2
	c2 = D2 / l^2 / h2

	u = @view x[1:n]
	v = @view x[n+1:2n]
	w = @view x[2n+1:3n]

	# Neumann boundary conditions
	f[1]   = v[1] - u[1]
	f[n]   = v[n] - u[n]

	f[n+1] = c1 * (	      - 1v[1] + v[2] ) + pp - (qq + 2) * v[1] + f1(v[1], w[1]) + u[1]
	f[2n]  = c1 * (v[n-1] - 1v[n]        ) + pp - (qq + 2) * v[n] + f1(v[n], w[n]) + u[n]

	f[n+1] = c2 * (       - 1w[1] + w[2])  + qq * w[1] - f1(v[1], w[1])
	f[end] = c2 * (w[n-1] - w[n]        )  + qq * w[n] - f1(v[n], w[n])

	for i=2:n-1
		f[i]    =                                  v[i] - u[i]
		f[n+i]  = c1 * (v[i-1] - 2v[i] + v[i+1]) + pp - (qq + 2) * v[i] + f1(v[i], w[i]) + u[i]
		f[2n+i] = c2 * (v[i-1] - 2v[i] + v[i+1]) + qq * w[i] - f1(v[i], w[i])
	end
	return f
end

function Laplacian1D(Nx, lx)
	hx = 2 * lx / Nx
	D2x = spdiagm(0 => -2 * ones(Nx), 1 => ones(Nx-1), -1 => ones(Nx-1)) / hx^2
	D2x[1, 1] = -1 / hx^2
	D2x[end, end] = -1 / hx^2

	D2xsp = sparse(D2x)
	return D2xsp
end

# Define the nonlinearites
p = 1.5;
q = 3.5;
f(u, v, w) = v - u
g(u, v, w) = p.*ones(length(v)) - (q+2)*v + v.^2 .* w + u
h(u, v, w) = q.*v-v.^2 .* w

# delta(n) = n == 0
# dϕ(u, κ, i, k) = (-1 + κ*exp(u[i]) / mean(exp.(u))) * delta(i-k) - κ*exp(u[i]) .* exp(u[k]) / (mean(exp.(u))*sum(exp.(u)))

function Fmit!(F, U, p)
	N = p.N
	u = U[1:N]
	v = U[N+1:2N]
	w = U[2N+1:3N]
	F[1:N] .= f(u, v, w)
	g1 = similar(v)
	mul!(g1, p.Δ, v)
	F[N+1:2N] .= p.diffcoef*g1 .+ g(u, v, w)
	h1 = similar(w)
	mul!(h1, p.Δ, w)
	F[2N+1:3N] .= 1*h1 .+ h(u, v, w)
	return F
end

# function JFmit(x,p)
# 	J = p.Δ
# 	dNonlin = zeros(length(x), length(x))
# 	for i in 1:length(x)
# 		for k in 1:length(x)
# 			dNonlin[i, k] = dϕ(x, p.κ, i, k)
# 		end
# 	end
# 	dNonlin = sparse(dNonlin)
# 	return p.diffcoef*J + dNonlin
# end


###########################################################################################################################
# Bifurcation Problem
# Parameters:
DiffCoef = 0.3;
Nx = 50; lx = 5;
Δ = Laplacian1D(Nx, lx);
par_mit = (diffcoef = DiffCoef, Δ = Δ, N = Nx);
par_bru = (pp = 1.5, qq = 3.5, D1 = DiffCoef, D2 = 1, l = 2*lx)
sol0 = vcat(par_bru.pp * ones(Nx), par_bru.pp * ones(Nx), par_bru.qq/par_bru.pp * ones(Nx))

# Define the L2 norm with weight.
# Choose the weighted norm in order to break symmetries and see all branches.
weight_norm = (lx .+ LinRange(-lx, lx, Nx)) |> vec
w_one = ones(Nx) |> vec
norm2(x) = norm(x .* w_one) / sqrt(length(x))
norm2_weighted(x) = norm(x .* weight_norm) / sqrt(length(x))
# initial guess f for newton
sol0 = ones(3*Nx) |> vec

nev_N = 6*Nx
int_param = [0.25, 0.35]

prob = BifurcationProblem(Fmit!, sol0, par_mit, (@optic _.diffcoef),;
  record_from_solution = (x, p; k...) -> (nrm = norm2(x[1:Int(end/3)]), nw = norm2_weighted(x[1:Int(end/3)]), n∞ = norminf(x[1:Int(end/3)]), sol=x),
  plot_solution = (x, p; k...) -> plot!(x[1:Int(end/3)] ; k...))


# eigensolver
eigls = EigArpack()

# options for Newton solver, we pass the eigen solver
opt_newton = BK.NewtonPar(tol = 1e-8, eigsolver = eigls, max_iterations = 20)

# options for continuation
opts_br = ContinuationPar(p_max = int_param[2], p_min = int_param[1],
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
	3,
	opts_br, bothside=true;
	kwargsC...
	)
p4 = plot(diagram, plotfold = true, markersize = 2, label = "", putspecialptlegend = false)
# ylims!(p4, 0, 5)  # Specify the y-axis limits for the plot

# If some branch in automatic bifurcation diagram is not computed, we can compute it by hand
# (1 = endpoint; 2,... = bifurcation points; length(diagram.γ.specialpoint) = endpoint):
br_missing = continuation(diagram.γ, 3,
		setproperties(opts_br; detect_bifurcation = 3, ds = 0.001, p_min = int_param[1], p_max = int_param[2]), bothside=true ;
		alg = PALC(),
		kwargsC...,
		nev = nev_N,
		)
plot!(p4, br_missing...)

# Plot the bifurcation diagram with weighted norm to see different branches	due to the symmetry breaking
p5 = plot(diagram; putspecialptlegend=false, markersize=2, vars = (:param, :nw), label = "",
		  title = "Plot with with weighted norm, #branches = $(size(diagram))")

######################################################################################################################
