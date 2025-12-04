using Revise, ForwardDiff
using BifurcationKit, LinearAlgebra, Plots, Statistics, SparseArrays
const BK = BifurcationKit

###############################################################################################################################

# Parameters
beta = [4.4, 1.2, 1.5, 4.0];
# nu2  = [0.0, 0.0002, 7.055];
nu2  = [0.0, 0.001, 0.5];
DiffCoef = 1.4;            # start diffusion for bifurcation parameter
N_species = 3;
int_param = [0.2, 2.5];    # Interval in which we consider bifurcation parameter
bif_param = 3;              # Number of diffusion which is used as bifurcation parameter

# Grid
L = 1.0;
N_fourier = 20;                 # number of Fourier modes
Nx = 100;               # number of collocation points for nonlinearities
x = range(0, stop=L, length=Nx);  # x_n = (n-1)/(N-1) for n=1, ..., N

# Laplace operator in Fourier space
kvec = 0:N_fourier;
laplace_factor = - (pi * kvec / L).^2;
# Cosinus-Matrix für N_fourier+1 Koeffizienten
C = [cos(k * π * x[n] / L) for n in 1:Nx, k in 0:N_fourier];  # C[n, k+1] = cos(k π x_n)

# Initial condition (Fourier coefficients)
# U0 = [42.59036855933033, 7981.33377368621, 16.2622428298969]; # for beta = [4.4, 1.2, 11.5, 4.8]
U0 = [4.4503021676176715, 87.1428332856512, 11.3534446020442]; # [for beta = [4.4, 1.2, 1.5, 4.0]
U0_real = hcat(U0[1] * ones(Nx), U0[2] * ones(Nx), U0[3] * ones(Nx));

# ---------------------------
# Helper functions
# ---------------------------
# Fourier -> real space reconstruction
function fourier_to_real(U_hat)
    Nx = size(C)[1]
    U_real = similar(U_hat, Nx, N_species)
    for s in 1:N_species
        U_real[:,s] = C * U_hat[:,s]
    end
    return U_real
end

function fourier_to_real_2(U_hat)
    N_real = size(C)[1]
    N_modes = size(U_hat,1)
    U_real = similar(U_hat, N_real, N_species)
    fill!(U_real, zero(eltype(U_hat)))   # WICHTIG: initialisieren
    for s in 1:N_species
        for k = 0:N_modes-1
            U_real[:,s] .+= U_hat[k+1,s] .* cos.(pi*k .* x ./ L)
        end
    end
    return U_real
end

function real_to_fourier(V)
    Nk = size(C,2)
    V_four = similar(V, Nk, N_species)
    # --- Least-Squares-Lösung ---
    for s in 1:N_species
        V_four[:,s] = C \ V[:,s]   # entspricht np.linalg.lstsq(C, V, rcond=None)[0]
    end
    return V_four
end

function test_fourier(V)
    V_hat = real_to_fourier(V)
    V_real = fourier_to_real(V_hat)
    V_real_2 = fourier_to_real_2(V_hat)

    V_rec = zeros(eltype(V), Nx, N_species)
    for s in 1:N_species
        for k in 0:N_fourier
            V_rec[:,s] .+= V_hat[k+1,s] .* cos.(pi * k .* x ./ L)
        end
        println("max Fehler = ", maximum(abs.(V[:,s] - V_real[:,s])))
        println("max Fehler_2 = ", maximum(abs.(V[:,s] - V_real_2[:,s])))
        println("max Fehler_rec = ", maximum(abs.(V[:,s] - V_rec[:,s])))
    end
end

# nonlinear function in real space
function F_nonl_real(U_real)
    Wl, Wd, C = U_real[:,1], U_real[:,2], U_real[:,3]
    f1 = beta[4]*Wd ./ ((1 .+ C) .* (1 .+ beta[2] .* Wl)) .- Wl
    f2 = beta[1] .* Wl .*Wl .- Wd
    f3 = Wd ./ (1 .+ beta[3] .* Wl) .- C
    return hcat(f1,f2,f3)
end

# Laplace applied to Fourier coefficients
# function apply_laplace(U_hat)
#     U_hat_new = similar(U_hat)
#     for s in 1:5
#         for k = 1:Nk
#             U_hat_new[k,s] = laplace_factor[k] * U_hat[k,s]
#         end
#     end
#     return U_hat_new
# end

# ---------------------------
# Nonlinearities in Fourier coefficients
# ---------------------------
function F_hat!(F_hat, U_hat, p)
    # 1. Transform to real space
    U_real = fourier_to_real(U_hat)
    
    # 2. Nonlinearities in real space
    Fnl_real = F_nonl_real(U_real)
    
    # 3. Project nonlinearities to Fourier space
    Fnl_hat = real_to_fourier(Fnl_real)
    
    # 4. Apply Laplacian in Fourier space
    for s in 1:N_species
        diffcoef = s==bif_param ? p.diffcoef : nu2[s]
        F_hat[:,s] .= laplace_factor .* U_hat[:,s] * diffcoef + Fnl_hat[:,s]
    end
    return F_hat
end

# ---------------------------
# Flattened nonlinearities for BifurcationKit
# ---------------------------
function F_flat!(F_flat, U_flat, p)
    # reshape flat vector to (Nk x 5) matrix
    Nk = size(C,2)
    U_hat = reshape(U_flat, Nk, N_species)
    F_hat = similar(U_hat)
    F_hat!(F_hat, U_hat, p)
    F_flat .= vec(F_hat)  # flatten
    return F_flat
end

# Initial solution in Fourier
U0_hat = real_to_fourier(U0_real)
# Initial solution flattened
sol0 = vec(U0_hat);

# Norms:
normC(x) = norm(x) / sqrt(length(x))
function nrm2_real(x_flat, ind_comp)
    Nk = size(C,2)
    x_unfl = reshape(x_flat, Nk, N_species)
    x_real = fourier_to_real(x_unfl)
    x_real_comp = x_real[:,ind_comp]
    return norm(x_real_comp)
end

function norminf_real(x_flat, ind_comp)
    Nk = size(C,2)
    x_unfl = reshape(x_flat, Nk, N_species)
    x_real = fourier_to_real(x_unfl)
    x_real_comp = x_real[:,ind_comp]
    return norminf(x_real_comp)
end

function norminf_all(x)
    return norminf_real(x,1)
end

# ---------------------------
# BifurcationKit setup
# ---------------------------
par_full = (diffcoef = DiffCoef,)
prob = BifurcationProblem(F_flat!, sol0, par_full, (@optic _.diffcoef);
                           record_from_solution = (x,p; k...) -> (nrmFirst=norm(x[1:Int(end/3)]), nrmReal=nrm2_real(x,1), nrm=norm(x), n∞ = norminf(x[1:Int(end/N_species)]), sol=x),
                           plot_solution = (x, p; k...) -> plot!(C * x[(bif_param-1)*Int(end/N_species)+1:bif_param*Int(end/N_species)] ; k...))
# options for Newton solver, we pass the eigen solver
# opt_newton = BK.NewtonPar(tol = 1e-10, max_iterations = 20);
opts_br = ContinuationPar(ds=1e-4, dsmax=5e-2, dsmin=1e-5,
                          p_min=int_param[1], p_max=int_param[2], nev=N_fourier, n_inversion=16,
                        #   newton_options=opt_newton,
                          detect_bifurcation=3, max_steps=400)

##############################################################################################################################
# Automatic Bifurcation diagram
diagram = @time bifurcationdiagram(prob, PALC(), 2, opts_br, bothside=true; verbosity=0, plot=true)
p1 = plot(diagram; markersize=2, title="Bifurcation diagram for nu_2 = $(nu2[2])", label="", vars = (:param, :nrmFirst))

# If some branch in automatic bifurcation diagram is not computed, we can compute it by hand
# (1 = endpoint; 2,... = bifurcation points; length(diagram.γ.specialpoint) = endpoint):
br_missing = continuation(diagram.child[1].γ, 2,
		opts_br,
		alg = PALC(),
		bothside=true;
		plot=true
		)
plot!(p1, br_missing, vars = (:param, :nrmFirst))

###############################################################################################################################
br = continuation(prob, PALC(), opts_br, bothside=true, normC = normC)

all_branches = Vector{Vector{Branch}}(undef, length(br.specialpoint)-2)
for (index, point) in pairs(br.specialpoint)
    @show index
    @show point.param
    if index == 1 || index == length(br.specialpoint)
        # Those are only endpoints of interval for κ and not bifurcations points, so skip them
        continue
    end
    branches = continuation(br, index,
        setproperties(opts_br), bothside=true ;
        alg = PALC(),
    )
    # if branches == Branch[]
    # 	@show "No branches found for index $(index), try smaller ds"
    # 	branches = continuation(br, index,
    # 	setproperties(opts_br; ds = 1e-5),
    # 	bothside=true ;
    # 	alg = PALC(),
    #     nev = 30,
    # 	)
    # end
    # normalize to always be Vector{Branch}
    if branches === nothing
        all_branches[index-1] = Branch[]  # empty vector
    elseif isa(branches, Branch)
        all_branches[index-1] = [branches]  # wrap single Branch in a vector
    elseif isa(branches, Vector{Branch})
        all_branches[index-1] = branches
    else
        error("Unexpected return type from continuation: $(typeof(branches))")
    end
end

# Plot main branch
p1 = plot(br; vars = (:param, :nrmFirst))

# Plot all secondary branches
for branch_vector in all_branches      # each element is Vector{Branch}
    for branch in branch_vector        # iterate individual Branch objects
        plot!(p1, branch; vars = (:param, :nrmFirst))
    end
end
display(p1)

