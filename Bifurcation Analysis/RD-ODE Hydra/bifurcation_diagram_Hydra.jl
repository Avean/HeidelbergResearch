using Revise, ForwardDiff
using BifurcationKit, LinearAlgebra, Plots, Statistics, SparseArrays
const BK = BifurcationKit

###############################################################################################################################

# Define the discretisation of our Problem
# Laplacian operator in 1D with Neumann boundary conditions

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

function Fmit!(F, U, p)
	N = p.N
	u = U[1:N]
	v = U[N+1:2N]
	w = U[2N+1:3N]
    f1 = similar(u)
	mul!(f1, p.Δ, u)
	F[1:N] .= 0.001*f1 .+ f(u, v, w)
	g1 = similar(v)
	mul!(g1, p.Δ, v)
	F[N+1:2N] .= p.diffcoef*g1 .+ g(u, v, w)
	h1 = similar(w)
	mul!(h1, p.Δ, w)
	F[2N+1:3N] .= 1.0*h1 .+ h(u, v, w)
	return F
end

# Jacobian of system:
# delta(n) = n == 0
# dϕ(u, κ, i, k) = (-1 + κ*exp(u[i]) / mean(exp.(u))) * delta(i-k) - κ*exp(u[i]) .* exp(u[k]) / (mean(exp.(u))*sum(exp.(u)))

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


############################################################################################################################
# Bifurcation Problem
# Parameters:
DiffCoef = 0.3;
Nx = 50; lx = 5;
Δ = Laplacian1D(Nx, lx);
par_mit = (diffcoef = DiffCoef, Δ = Δ, N = Nx);
sol0 = vcat(par_bru.pp * ones(Nx), par_bru.pp * ones(Nx), par_bru.qq/par_bru.pp * ones(Nx))

# Define the L2 norm with weight.
# Choose the weighted norm in order to break symmetries and see all branches.
weight_norm = (lx .+ LinRange(-lx, lx, Nx)) |> vec
w_one = ones(Nx) |> vec
norm2(x) = norm(x .* w_one) / sqrt(length(x))
norm2_weighted(x) = norm(x .* weight_norm) / sqrt(length(x))

nev_N = 50 #6*Nx
int_param = [0.25, 0.35]

prob = BifurcationProblem(Fmit!, sol0, par_mit, (@optic _.diffcoef),;
  record_from_solution = (x, p; k...) -> (nrm = norm2(x[1:Int(end/3)]), nw = norm2_weighted(x[1:Int(end/3)]), n∞ = norminf(x[1:Int(end/3)]), sol=x),
  plot_solution = (x, p; k...) -> plot!(x[1:Int(end/3)] ; k...))


# eigensolver
eigls = EigArpack()

# options for Newton solver, we pass the eigen solver
opt_newton = BK.NewtonPar(tol = 1e-8, verbose = true, eigsolver = eigls, max_iterations = 20)

# options for continuation
opts_br = ContinuationPar(p_min = int_param[1], p_max = int_param[2],
	# for a good looking curve
	dsmin = 0.001, dsmax = 0.01, ds = 0.001,
	# detect codim 1 bifurcations
	detect_bifurcation = 3,
    # number of eigenvalues to compute
	nev = nev_N,
    plot_every_step = 10,
	newton_options = (@set opt_newton.verbose = false),
    max_steps = 251,
    tol_stability = 1e-6,
	n_inversion = 6,
    # Optional: bisection options for locating bifurcations
    dsmin_bisection = 1e-7, max_bisection_steps = 25, tol_bisection_eigenvalue = 1e-19
    )

############################################################################################################################
# Calculating branches one by one.

br = continuation(prob, PALC(), opts_br, bothside=true, normC = norminf)

all_branches = Array{Branch}(undef, length(br.specialpoint)-2)
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
		nev = nev_N,
	)
	if branches == Branch[]
		@show "No branches found for index $(index), try smaller ds"
		branches = continuation(br, index,
		setproperties(opts_br; ds = 0.001),
		bothside=true ;
		alg = PALC(),
		nev = nev_N,
		)
	end
	all_branches[index-1] = branches
end

p1 = plot(br)
for branch in all_branches
	p1 = plot!(branch)
end
display(p1)

# Use weighted norm for plotting to see different branches due to the symmetry breaking
p2 = plot(br; vars = (:param, :nw), title = "Weighted norm")
for branch in all_branches
	p2 = plot!(branch; vars = (:param, :nw))
end
plot(p1, p2)

######################################################################################################################
# One can also use bifurcationdiagram() to compute the bifurcation diagram automatically and then by hand the missing branches.
# automatic bifurcation diagram computation
diagram = @time bifurcationdiagram(prob, PALC(),
	# very important parameter. This specifies the maximum amount of recursion
	# when computing the bifurcation diagram. It means we allow computing branches of branches
	# at most in the present case.
	3,
	opts_br, bothside=true;
    verbosity = 0, plot = true,
    # callback_newton = cb,
	usedeflation = true,
	# finalise_solution = finSol,
	normC = norminf
	)
p4 = plot(diagram; plotfold = false, putspecialptlegend=false, markersize = 2, title = "#branches = $(size(diagram))", label="")
# ylims!(p4, 0, 5)  # Specify the y-axis limits for the plot

# If some branch in automatic bifurcation diagram is not computed, we can compute it by hand
# (1 = endpoint; 2,... = bifurcation points; length(diagram.γ.specialpoint) = endpoint):
br_missing = continuation(diagram.γ, 3,
		setproperties(opts_br; detect_bifurcation = 3, ds = 0.001, p_min = int_param[1], p_max = int_param[2]), bothside=true ;
		alg = PALC(),
		normC = norminf,
		nev = nev_N,
		)
plot!(p4, br_missing)

# Continue already computed bifurcation diagram:
bifurcationdiagram!(prob,
	# this improves the first branch on the red? curve. Note that
	# for symmetry reasons, the first bifurcation point
	# has ? branches
	get_branch(diagram, (1,)), 5, opts_br;
	verbosity = 0, plot = true,
	# callback_newton = cb,
	# finalise_solution = finSol,
	usedeflation = true,
	normC = norminf)

# Plot the bifurcation diagram with weighted norm to see different branches	due to the symmetry breaking
p5 = plot(diagram; putspecialptlegend=false, markersize=2, vars = (:param, :nw), label="",
		  title = "Plot with with weighted norm, #branches = $(size(diagram))")
# Plot the bifurcation diagram with supremum norm
p6 = plot(diagram; putspecialptlegend=false, markersize=2, vars = (:param, :n∞), label="",
		  title = "Plot with with supremum norm, #branches = $(size(diagram))")

######################################################################################################################
