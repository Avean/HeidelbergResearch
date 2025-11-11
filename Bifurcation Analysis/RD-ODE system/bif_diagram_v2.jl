using Revise
using BifurcationKit, Plots, SparseArrays
const BK = BifurcationKit

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

function Jbru_sp(x, p)
	(;pp, qq, D1, D2, l) = p
	# compute the Jacobian using a sparse representation
	n = div(length(x), 3)
	h = 1.0 / n; h2 = h*h

	c1 = D1 / p.l^2 / h2
	c2 = D2 / p.l^2 / h2

	u = @view x[1:n]
	v = @view x[n+1:2n]
	w = @view x[2n+1:3n]

	diag   = zeros(eltype(x), 3n)
	diagp1 = zeros(eltype(x), 3n-1)
	diagm1 = zeros(eltype(x), 3n-1)

	diagpn = zeros(eltype(x), n)
	diagmn = zeros(eltype(x), n)

	@. diagmn = β - 2 * u * v
	@. diagm1[n:2n-1] = c1
	@. diagm1[2n+1:end] = c2

	@. diag[1:n]    = 1
	@. diag[n+1:2n] = -2c1 - (qq + 2) + 2 * v * w
	@. diag[2n+1:3n] = -2c2 + qq - 2 * v * w

	@. diagp1[n:2n-1] = c1
	@. diagp1[2n+1:end] = c2

	@. diagpn = u * u
	return spdiagm(0 => diag, 1 => diagp1, -1 => diagm1, n => diagpn, -n => diagmn)
end

n = Int(100*3)

# parameters of the Brusselator model and guess for the stationary solution
par_bru = (pp = 1.5, qq = 3.5, D1 = 0.3, D2 = 1, l = 10)
sol0 = vcat(par_bru.pp * ones(Int(n/3)), par_bru.pp * ones(Int(n/3)), par_bru.qq/par_bru.pp * ones(Int(n/3)))
int_param = [0.25, 0.35]

# bifurcation problem
probBif = BK.BifurcationProblem(Fbru!, sol0, par_bru, (@optic _.D1);
#   J = Jbru_sp,
  plot_solution = (x, p; kwargs...) -> plot!(x[1:Int(end/3)] ; kwargs...),
  record_from_solution = (x, p; k...) -> x[div(n,3)])

eigls = EigArpack(1.1, :LM)

opt_newton = NewtonPar(eigsolver = eigls, tol = 1e-9)
opts_br_eq = ContinuationPar(dsmin = 0.001, dsmax = 0.01, ds = 0.001,
	p_max = int_param[2], p_min = int_param[1], nev = 21,
	newton_options = opt_newton, max_steps = 1000,
	# specific options for precise localization of Hopf points
	n_inversion = 6)

kwargsC = (verbosity = 0, plot = true, normC = norminf)

br = continuation(probBif, PALC(), opts_br_eq, normC = norminf)

scene = plot(br)

# automatic bifurcation diagram generation
diagram = @time bifurcationdiagram(probBif, PALC(),
	# very important parameter. This specifies the maximum amount of recursion
	# when computing the bifurcation diagram. It means we allow computing branches of branches
	# at most in the present case.
	2,
	opts_br_eq, bothside=true;
	kwargsC...,
	)
p4 = plot(diagram, plotfold = true, markersize = 2, label = "", putspecialptlegend = false)

# automatic branch switching from Hopf point
opt_po = NewtonPar(tol = 1e-10, verbose = true, max_iterations = 15)
opts_po_cont = ContinuationPar(dsmin = 0.001,
		dsmax = 0.04, ds = 0.01,
		p_max = 2.2,
		max_steps = 30,
		newton_options = opt_po,
		plot_every_step = 1,
		nev = 11,
		tol_stability = 1e-6,
		)

# number of time slices for the periodic orbit
M = 51
probFD = PeriodicOrbitTrapProblem(M = M;
  # specific method for solving the linear system for newton
  # of periodic orbits with trapeze method.
  # We could use the default one :FullLU (slower here)
  jacobian = :BorderedSparseInplace)

br_po = continuation(
	# arguments for branch switching from the first
	# Hopf bifurcation point
	br, 1,
	# arguments for continuation
	opts_po_cont, probFD;
	# regular options for continuation
	verbosity = 3, plot = true,
	plot_solution = (x, p; kwargs...) -> heatmap!(reshape(x[1:end-1], 2*n, M)'; ylabel="time", color=:viridis, kwargs...),
	normC = norminf)

Scene = title!("")

br_po2 = continuation(
	# arguments for branch switching
	br_po, 1,
	# arguments for continuation
	opts_po_cont;
	ampfactor = 1., δp = 0.01,
	verbosity = 3, plot = true,
	plot_solution = (x, p; kwargs...) -> heatmap!(reshape(x[1:end-1], 2*n, M)'; ylabel="time", color=:viridis, kwargs...),
	normC = norminf)

