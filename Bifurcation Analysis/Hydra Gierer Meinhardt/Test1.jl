using Revise
using BifurcationKit, Plots, SparseArrays
const BK = BifurcationKit

f1(u, v) = u * u * v

function Fbru!(f, x, p, t = 0)
	(;α, β, D1, D2, l) = p
	n = div(length(x), 2)
	h2 = 1.0 / n^2
	c1 = D1 / l^2 / h2
	c2 = D2 / l^2 / h2

	u = @view x[1:n]
	v = @view x[n+1:2n]

	# Dirichlet boundary conditions
	f[1]   = c1 * (α	  - 2u[1] + u[2] ) + α - (β + 1) * u[1] + f1(u[1], v[1])
	f[end] = c2 * (v[n-1] - 2v[n] + β / α)			 + β * u[n] - f1(u[n], v[n])

	f[n]   = c1 * (u[n-1] - 2u[n] +  α   ) + α - (β + 1) * u[n] + f1(u[n], v[n])
	f[n+1] = c2 * (β / α  - 2v[1] + v[2])			 + β * u[1] - f1(u[1], v[1])

	for i=2:n-1
		  f[i] = c1 * (u[i-1] - 2u[i] + u[i+1]) + α - (β + 1) * u[i] + f1(u[i], v[i])
		f[n+i] = c2 * (v[i-1] - 2v[i] + v[i+1])			  + β * u[i] - f1(u[i], v[i])
	end
	return f
end

function Jbru_sp(x, p)
	(;α, β, D1, D2, l) = p
	# compute the Jacobian using a sparse representation
	n = div(length(x), 2)
	h = 1.0 / n; h2 = h*h

	c1 = D1 / p.l^2 / h2
	c2 = D2 / p.l^2 / h2

	u = @view x[1:n]
	v = @view x[n+1:2n]

	diag   = zeros(eltype(x), 2n)
	diagp1 = zeros(eltype(x), 2n-1)
	diagm1 = zeros(eltype(x), 2n-1)

	diagpn = zeros(eltype(x), n)
	diagmn = zeros(eltype(x), n)

	@. diagmn = β - 2 * u * v
	@. diagm1[1:n-1] = c1
	@. diagm1[n+1:end] = c2

	@. diag[1:n]    = -2c1 - (β + 1) + 2 * u * v
	@. diag[n+1:2n] = -2c2 - u * u

	@. diagp1[1:n-1] = c1
	@. diagp1[n+1:end] = c2

	@. diagpn = u * u
	return spdiagm(0 => diag, 1 => diagp1, -1 => diagm1, n => diagpn, -n => diagmn)
end

n = 300

# parameters of the Brusselator model and guess for the stationary solution
par_bru = (α = 2., β = 5.45, D1 = 0.008, D2 = 0.004, l = 0.3)
sol0 = vcat(par_bru.α * ones(n), par_bru.β/par_bru.α * ones(n))

# bifurcation problem
probBif = BK.BifurcationProblem(Fbru!, sol0, par_bru, (@optic _.l);
  J = Jbru_sp,
  plot_solution = (x, p; kwargs...) -> (plotsol(x; label="", kwargs... )),
  record_from_solution = (x, p; k...) -> x[div(n,2)])

  eigls = EigArpack(1.1, :LM)

  opt_newton = NewtonPar(eigsolver = eigls, tol = 1e-9)
opts_br_eq = ContinuationPar(dsmin = 0.001, dsmax = 0.01, ds = 0.001,
	p_max = 1.9, nev = 21,
	newton_options = opt_newton, max_steps = 1000,
	# specific options for precise localization of Hopf points
	n_inversion = 6)

br = continuation(probBif, PALC(), opts_br_eq, normC = norminf)

scene = plot(br)

hopfpt = get_normal_form(br, 1)

# index of the Hopf point in br.specialpoint
ind_hopf = 2

# newton iterations to compute the Hopf point
hopfpoint = newton(br, ind_hopf; normN = norminf)
BK.converged(hopfpoint) && printstyled(color=:red, "--> We found a Hopf Point at l = ", hopfpoint.u.p[1], ", ω = ", hopfpoint.u.p[2], ", from l = ", br.specialpoint[ind_hopf].param, "\n")



optcdim2 = ContinuationPar(dsmin = 0.001, dsmax = 0.05, ds= 0.01, p_max = 6.5, p_min = 0.0, newton_options = opt_newton, detect_bifurcation = 0)

br_hopf = continuation(br, ind_hopf, (@optic _.β),
	optcdim2, verbosity = 2,
	# detection of codim 2 bifurcations with bisection
	detect_codim2_bifurcation = 2,
	# we update the Hopf problem at every continuation step
	update_minaug_every_step = 1,
	jacobian_ma = BK.MinAug(), # specific to large dimensions
	normC = norminf,
	)

scene = plot(br_hopf)


###

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
  # We could use the default one FullLU() (slower here)
  jacobian = BK.BorderedSparseInplace())

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


###

br_po2 = continuation(
	# arguments for branch switching
	br_po, 1,
	# arguments for continuation
	opts_po_cont;
	ampfactor = 1., δp = 0.01,
	verbosity = 3, plot = true,
	plot_solution = (x, p; kwargs...) -> heatmap!(reshape(x[1:end-1], 2*n, M)'; ylabel="time", color=:viridis, kwargs...),
	normC = norminf)