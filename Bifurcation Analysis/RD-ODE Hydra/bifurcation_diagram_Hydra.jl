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
nu = [0.0 3.8154e-05 0.4433 6.0713e-08 0.0004];
beta = [1.0629 540.4003 1.1596 11.5964 11.5964 4.8254];
F_nonl(Wl, A, Wd, C, S) = [beta[6]*S ./ ((1.0 .+ A).*(1.0 .+ C).*(1.0 .+ beta[3]*Wl)) - Wl,
                           beta[1] ./ (1.0 .+ beta[4]*Wl) - A,
                           beta[2]*Wl.*S - Wd,
                           Wd ./ (1.0 .+ beta[5]*Wl) - C,
                           Wl - S
                          ]

# Non-trivial constant steady state:
U0 = [0.08167547924306925, 0.54587711524379, 3.6049476660047937, 1.8514050545898464, 0.08167547924306925]

function Fmit!(F, U, p)
	N = p.N
	Wl = U[1:N]
	A = U[N+1:2N]
	Wd = U[2N+1:3N]
    C = U[3N+1:4N]
    S = U[4N+1:5N]

    f1 = similar(Wl)
	mul!(f1, p.Δ, Wl)
	F[1:N] .= nu[1]*f1 .+ F_nonl(Wl, A, Wd, C, S)[1]

	f2 = similar(A)
	mul!(f2, p.Δ, A)
	F[N+1:2N] .= nu[2]*f2 .+ F_nonl(Wl, A, Wd, C, S)[2]
    
    f3 = similar(Wd)
	mul!(f3, p.Δ, Wd)
	F[2N+1:3N] .= p.diffcoef*f3 .+ F_nonl(Wl, A, Wd, C, S)[3]

    f4 = similar(C)
	mul!(f4, p.Δ, C)
	F[3N+1:4N] .= nu[4]*f4 .+ F_nonl(Wl, A, Wd, C, S)[4]

    f5 = similar(S)
	mul!(f5, p.Δ, S)
	F[4N+1:5N] .= nu[5]*f5 .+ F_nonl(Wl, A, Wd, C, S)[5]

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
DiffCoef = 0.009;
Nx = 50; lx = 0.5;
Δ = Laplacian1D(Nx, lx);
par_mit = (diffcoef = DiffCoef, Δ = Δ, N = Nx);
sol0 = vcat(U0[1] * ones(Nx), U0[2] * ones(Nx), U0[3] * ones(Nx), U0[4] * ones(Nx), U0[5] * ones(Nx))

# Define the L2 norm with weight.
# Choose the weighted norm in order to break symmetries and see all branches.
weight_norm = (lx .+ LinRange(-lx, lx, Nx)) |> vec
w_one = ones(Nx) |> vec
norm2(x) = norm(x .* w_one) / sqrt(length(x))
norm2_weighted(x) = norm(x .* weight_norm) / sqrt(length(x))

nev_N = 5*Nx
int_param = [0.001, 0.5]

prob = BifurcationProblem(Fmit!, sol0, par_mit, (@optic _.diffcoef),;
  record_from_solution = (x, p; k...) -> (nrm = norm2(x[1:Int(end/5)]), nw = norm2_weighted(x[1:Int(end/5)]), n∞ = norminf(x[1:Int(end/5)]), sol=x),
  plot_solution = (x, p; k...) -> plot!(x[1:Int(end/5)] ; k...))


# eigensolver
eigls = EigArpack()

# options for Newton solver, we pass the eigen solver
opt_newton = BK.NewtonPar(tol = 1e-8, verbose = true, eigsolver = eigls, max_iterations = 20)

# options for continuation
opts_br = ContinuationPar(p_min = int_param[1], p_max = int_param[2],
	# for a good looking curve
	dsmin = 0.0001, dsmax = 0.05, ds = 0.005,
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
	2,
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
