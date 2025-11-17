using Revise, ForwardDiff
using BifurcationKit, LinearAlgebra, Plots, Statistics, SparseArrays
const BK = BifurcationKit

includet("../../Evolution Solvers/Models/GiereMeinhardt/GiereMeinhardtModules.jl")

includet("../../Evolution Solvers/FillFunctions.jl")
includet("../../Evolution Solvers/DiffusionMatrices.jl")
###############################################################################################################################

# Define the discretisation of our Problem
# Laplacian operator in 1D with Neumann boundary conditions

using .LaplaceDiscretisation
using .Struktury


# Define the nonlinearites
nu = [0.0, 3.8154e-05, 0.4433, 6.0713e-08, 0.0004];
beta = [1.0629 540.4003 1.1596 11.5964 11.5964 4.8254];
F_nonl(Wl, A, Wd, C, S) = [beta[6]*S ./ ((1.0 .+ A).*(1.0 .+ C).*(1.0 .+ beta[3]*Wl)) - Wl,
                           beta[1] ./ (1.0 .+ beta[4]*Wl) - A,
                           beta[2]*Wl.*S - Wd,
                           Wd ./ (1.0 .+ beta[5]*Wl) - C,
                           Wl - S
                          ]

# Non-trivial constant steady state:
U0 = [0.08167547924306925, 0.54587711524379, 3.6049476660047937, 1.8514050545898464, 0.08167547924306925];

Nodes = 15;
FullMatrix = SC.N.Full[1:Nodes,:] .* sqrt(SimParam.N-1);
TrunMatrix = SC.N.Trunc[1:Nodes,:] .* sqrt(SimParam.N-1);
Eigen = Eig.C[1:Nodes];

function FGMFourier!(f, x, DiffCoefBif, t=0)
    
    Dom = range(0, 1, SimParam.N)
    dx = 1/(SimParam.N-1);
    
    n = div(length(x), 5)

    WntLOCCoeff     = @view x[1:n]
    DkkACoeff       = @view x[n+1:2*n]
    WntDiffCoeff    = @view x[2*n+1:3*n]
    DkkCCoeff       = @view x[3*n+1:4*n]
    SDCoeff         = @view x[4*n+1:5*n]

    WntLOC  = FullMatrix' * WntLOCCoeff
    DkkA    = FullMatrix' * DkkACoeff
    WntDiff = FullMatrix' * WntDiffCoeff
    DkkC    = FullMatrix' * DkkCCoeff
    SD      = FullMatrix' * SDCoeff

    WntLOCaCoeff     = TrunMatrix * (par_GM.β6 .* SD ./ (1 .+ DkkA) ./ (1.0 .+ DkkC) ./ (1.0 .+ par_GM.β3*WntLOC) .- WntLOC) ./ (SimParam.N - 1)
    DkkAaCoeff       = TrunMatrix * (par_GM.β1 ./ (1 .+ par_GM.β4.*WntLOC) .- DkkA) ./ (SimParam.N - 1)
    WntDiffaCoeff    = TrunMatrix * (par_GM.β2 .* WntLOC.*SD .- WntDiff) ./ (SimParam.N - 1)
    DkkCaCoeff       = TrunMatrix * (WntDiff ./ (1.0 .+ par_GM.β5 .* WntLOC) .- DkkC) ./ (SimParam.N - 1)
    SDaCoeff         = TrunMatrix * (WntLOC .- SD) ./ (SimParam.N - 1) 

    for i in 1:n
        f[i]        = - par_GM.ν1      * Eigen[i] * WntLOCCoeff[i]   + WntLOCaCoeff[i]
        f[n+i]      = - par_GM.ν2      * Eigen[i] * DkkACoeff[i]     + DkkAaCoeff[i]
        f[2*n+i]    = - DiffCoefBif.ν3 * Eigen[i] * WntDiffCoeff[i]  + WntDiffaCoeff[i]
        f[3*n+i]    = - par_GM.ν4      * Eigen[i] * DkkCCoeff[i]     + DkkCaCoeff[i]
        f[4*n+i]    = - par_GM.ν5      * Eigen[i] * SDCoeff[i]      + SDaCoeff[i]
    end

    return f
end


ν = [0.0, 3.8154e-05, 0.4433, 6.0713e-08, 0.0004];
β = [1.0629, 540.4003, 1.1596, 11.5964, 11.5964, 4.8254];

ν[3] = 0.03
par_GM = (
    ν1 = ν[1],
    ν2 = ν[2],
    ν3 = ν[3],
    ν4 = ν[4],
    ν5 = ν[5],
    β1 = β[1],
    β2 = β[2],
    β3 = β[3],
    β4 = β[4],
    β5 = β[5],
    β6 = β[6]
)

DiffCoefBif = (ν3 =  0.22, Da = 0.01)

############################################################################################################################
# Bifurcation Problem

sol0 = zeros(5*Nodes)
sol0[1] = U0[1]
sol0[11] = U0[2]
sol0[21] = U0[3]
sol0[31] = U0[4]
sol0[41] = U0[5]


Nx = 15; lx = 0.5;
# Define the L2 norm with weight.
# Choose the weighted norm in order to break symmetries and see all branches.
weight_norm = (lx .+ LinRange(-lx, lx, Nx)) |> vec
w_one = ones(Nx) |> vec
norm2(x) = norm(x .* w_one) / sqrt(length(x));
norm2_weighted(x) = norm(x .* weight_norm) / sqrt(length(x));

# unweighted L2 over all species
normC(x) = norm(x) / sqrt(length(x))
# # weighted L2 emphasizing Wl
# weight = ones(5*Nx)
# weight[1:Nx] .= 1.0  # Wl
# weight[Nx+1:5Nx] .= 0.5  # other species
# normC(x) = norm(x .* weight) / sqrt(length(x))



prob = BifurcationProblem(FGMFourier!, sol0, DiffCoefBif, (@optic _.ν3),;
	# J = jacobian,
	record_from_solution = (x, p; k...) -> (nrm = norm2(x[1:Int(end/5)]), nw = norm2_weighted(x[1:Int(end/5)]), n∞ = norminf(x[1:Int(end/5)]), sol=x),
	plot_solution = (x, p; k...) -> plot!(x[1:Int(end/5)] ; k...))


int_param = [0.02, 0.3]  # interval of continuation for the diffusion coefficient

# Beispielwerte, anpassbar
nev_N = 15; # number of eigenvalues to compute
eig_ncv = min(5*nev_N, 5*Nx*5);   # aber ≤ total dim; choose sensible cap
eig_ncv = min(eig_ncv, 5*Nx);     # ensure <= Ntot
# eig_ncv = 40;
eig_tol = 1e-6;
eig_maxiter = 2000;

# eigensolver
eigls = EigArpack(ncv = eig_ncv, tol = eig_tol, maxiter = eig_maxiter); #ncv = eig_ncv, tol = eig_tol, maxiter = eig_maxiter
# options for Newton solver, we pass the eigen solver
opt_newton = BK.NewtonPar(tol = 1e-8, verbose = true, eigsolver = eigls, max_iterations = 20);



# options for continuation
opts_br = ContinuationPar(p_min = int_param[1], p_max = int_param[2],
	# for a good looking curve
	dsmin = 1e-4, dsmax = 1e-2, ds = 1e-3,
	# detect codim 1 bifurcations
	detect_bifurcation = 3,  # reduziert Kosten
    # number of eigenvalues to compute
	nev = nev_N,
    plot_every_step = 20, 
	newton_options = (@set opt_newton.verbose = false),
    # max_steps = 200,
    # tol_stability = 1e-6,
	# n_inversion = 6,
    # Optional: bisection options for locating bifurcations
    # dsmin_bisection = 1e-7, max_bisection_steps = 25, tol_bisection_eigenvalue = 1e-8
    );

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
	usedeflation = false,
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