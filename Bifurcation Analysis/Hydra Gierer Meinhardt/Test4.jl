using Revise, ForwardDiff
using BifurcationKit, LinearAlgebra, Plots, Statistics, SparseArrays
const BK = BifurcationKit

###############################################################################################################################

# Define the discretisation of our Problem
# Laplacian operator in 1D with Neumann boundary conditions

function Laplacian1D(Nx, lx)
	hx = 2*lx / (Nx-1) # hx = 2 * lx / Nx
	D2x = spdiagm(0 => -2 * ones(Nx), 1 => ones(Nx-1), -1 => ones(Nx-1)) / hx^2
	D2x[1, 1] = -1 / hx^2
	D2x[end, end] = -1 / hx^2

	D2xsp = sparse(D2x)
	return D2xsp
end

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


function FGMeinhardt!(f, x, DiffCoefBif, t = 0)

    n = div(length(x), 5)
    dx = 1/(n-1);

    WntLOC =  @view x[1:n]
    DkkA =    @view x[n+1:2*n]
    WntDiff = @view x[2*n+1:3*n]
    DkkC =    @view x[3*n+1:4*n]
    SD =      @view x[4*n+1:5*n]
	
    WntLOCa = par_GM.β6 .* SD ./ (1 .+ DkkA) ./ (1.0 .+ DkkC) ./ (1.0 .+ par_GM.β3*WntLOC) .- WntLOC
    DkkAa = par_GM.β1 ./ (1 .+ par_GM.β4.*WntLOC) .- DkkA
    WntDiffa = par_GM.β2 .* WntLOC.*SD .- WntDiff
    DkkCa = WntDiff ./ (1.0 .+ par_GM.β5 .* WntLOC) .- DkkC
    SDa = WntLOC .- SD

    WntLOCb =   similar(WntLOC)
    DkkAb =     similar(DkkA)
    WntDiffb =  DiffCoefBif.ν3 * similar(WntDiff)
    DkkCb =     similar(DkkC)
    SDb =       similar(SD)
    
    # println("To typ Diff = ", typeof(WntDiff))
    # print("To typ Diffb = ", typeof(WntDiffb))
    # # sleep(2)
    # println()
    # println()
    # println()

    # WntLOCb1  = par_GM.ν1 * Δ * WntLOC 
    # DkkAb1    = par_GM.ν2 * Δ * DkkA 
    # WntDiffb1 = DiffCoefBif.ν3 * Δ * WntDiff
    # DkkCb1    = par_GM.ν4 * Δ * DkkC
    # SDb1      = par_GM.ν5 * Δ * SD                 
    
    Lap2(V,i) = V[i-1] + V[i+1] - 2*V[i]

    for i in 2:n-1
        WntLOCb[i]  = par_GM.ν1 * Lap2(WntLOC,i)/dx^2
        DkkAb[i]    = par_GM.ν2 * Lap2(DkkA,i)/dx^2
        WntDiffb[i] = DiffCoefBif.ν3 * Lap2(WntDiff,i)/dx^2;
        DkkCb[i]    = par_GM.ν4 * Lap2(DkkC,i)/dx^2;
        SDb[i]      = par_GM.ν5 * Lap2(SD,i)/dx^2;
    end

    # Boundary conditions for all 5 variables
    WntLOCb[1] = par_GM.ν1 * (WntLOC[2] - WntLOC[1]) / dx^2
    WntLOCb[end] = par_GM.ν1 * (WntLOC[end-1] - WntLOC[end]) / dx^2

    DkkAb[1] = par_GM.ν2 * (DkkA[2] - DkkA[1]) / dx^2
    DkkAb[end] = par_GM.ν2 *  (DkkA[end-1] - DkkA[end]) / dx^2

    WntDiffb[1] = DiffCoefBif.ν3 * (WntDiff[2] - WntDiff[1]) / dx^2
    WntDiffb[end] = DiffCoefBif.ν3 * (WntDiff[end-1] - WntDiff[end]) / dx^2

    DkkCb[1] = par_GM.ν4 * (DkkC[2] - DkkC[1]) / dx^2
    DkkCb[end] = par_GM.ν4 * (DkkC[end-1] - DkkC[end]) / dx^2

    SDb[1] = par_GM.ν5 * (SD[2] - SD[1]) / dx^2
    SDb[end] = par_GM.ν5 * (SD[end-1] - SD[end]) / dx^2

    # println(ForwardDiff.value(norm(WntLOCb - WntLOCb1)))
    # println(ForwardDiff.value(norm(DkkAb - DkkAb1)))
    # println(ForwardDiff.value(norm(WntDiffb - WntDiffb1)))
    # println(ForwardDiff.value(norm(DkkCb - DkkCb1)))
    # println(ForwardDiff.value(norm(SDb - SDb1)))

    # println()
    # println()
    # println()


    # Assemble f for all 5 variables
    for i in 1:n
        f[i] = WntLOCb[i] + WntLOCa[i]
        f[n+i] = DkkAb[i] + DkkAa[i]
        f[2*n+i] = WntDiffb[i] + WntDiffa[i]
        f[3*n+i] = DkkCb[i] + DkkCa[i]
        f[4*n+i] = SDb[i] + SDa[i]
    end

    # f[1:n]       = WntLOCb + WntLOCa
    # f[n+1:2*n]   = DkkAb + DkkAa
    # f[2*n+1:3*n] = WntDiffb + WntDiffa
    # f[3*n+1:4*n] = DkkCb + DkkCa
    # f[4*n+1:5*n] = SDb + SDa
	
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

DiffCoefBif = (ν3 =  0.03, Da = 0.01)

# assume:
#   nu  :: Vector{Float64} length 5
#   beta:: Vector{Float64} length 6
#   p.Δ :: SparseMatrixCSC  (Nx×Nx Laplacian)
#   p.N :: Int
#   p.diffcoef :: Float64   (diffusion coefficient of S)


############################################################################################################################
# Bifurcation Problem
# Parameters:
DiffCoef = 0.03;
Nx = 15; lx = 0.5;
Δ = Laplacian1D(Nx, lx);
par_mit = (diffcoef = DiffCoef, Δ = Δ, N = Nx);
sol0 = vcat(U0[1] * ones(Nx), U0[2] * ones(Nx), U0[3] * ones(Nx), U0[4] * ones(Nx), U0[5] * ones(Nx))

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


prob = BifurcationProblem(Fmit!, sol0, par_mit, (@optic _.diffcoef),;
	# J = jacobian,
	record_from_solution = (x, p; k...) -> (nrm = norm2(x[1:Int(end/5)]), nw = norm2_weighted(x[1:Int(end/5)]), n∞ = norminf(x[1:Int(end/5)]), sol=x),
	plot_solution = (x, p; k...) -> plot!(x[1:Int(end/5)] ; k...))

prob = BifurcationProblem(FGMeinhardt!, sol0, par_GM, (@optic _.ν3),;
	# J = jacobian,
	record_from_solution = (x, p; k...) -> (nrm = norm2(x[1:Int(end/5)]), nw = norm2_weighted(x[1:Int(end/5)]), n∞ = norminf(x[1:Int(end/5)]), sol=x),
	plot_solution = (x, p; k...) -> plot!(x[1:Int(end/5)] ; k...))


int_param = [0.018, 0.3]  # interval of continuation for the diffusion coefficient

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
    max_steps = 200,
    tol_stability = 1e-6,
	n_inversion = 6,
    # Optional: bisection options for locating bifurcations
    dsmin_bisection = 1e-7, max_bisection_steps = 25, tol_bisection_eigenvalue = 1e-8
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
    # callback_newton = cb,
	usedeflation = false,
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