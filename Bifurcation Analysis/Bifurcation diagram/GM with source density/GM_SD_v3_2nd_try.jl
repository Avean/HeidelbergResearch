using BifurcationKit, LinearAlgebra, Plots
const BK = BifurcationKit

# ============================================================
# PARAMETERS (keep Float64 ONLY)
# ============================================================

param_a = 2.0
param_b = 1.0
param_muu = 0.5
param_muuv = 1.0

nu2 = [0.005, 0.05, 1.0]
bif_param = 2
N_species = 3

L = 1.0
N_fourier = 20
Nx = 100

x = range(0, stop=L, length=Nx)
kvec = 0:N_fourier

laplace_factor = - (pi * kvec / L).^2
C = [cos(k * π * x[n] / L) for n in 1:Nx, k in 0:N_fourier]

# ============================================================
# TRANSFORMS (AD-safe)
# ============================================================

function fourier_to_real(U_hat)
    T = eltype(U_hat)
    C_T = T.(C)

    Nx = size(C,1)
    U_real = Array{T}(undef, Nx, N_species)

    for s in 1:N_species
        @views U_real[:,s] = C_T * U_hat[:,s]
    end
    return U_real
end

function real_to_fourier(V)
    T = eltype(V)
    C_T = T.(C)

    Nk = size(C,2)
    V_hat = Array{T}(undef, Nk, N_species)

    for s in 1:N_species
        @views V_hat[:,s] = C_T \ V[:,s]
    end
    return V_hat
end

# ============================================================
# NONLINEARITY (AD-safe)
# ============================================================

function F_nonl_real(U)
    T = eltype(U)

    S = U[:,1]
    Uu = U[:,2]
    Vv = U[:,3]

    f1 = Uu .- S
    f2 = T(param_a) .* Uu .* Uu .* S ./ (one(T) .+ Vv) .- T(param_muu) .* Uu
    f3 = T(param_b) .* Uu .* Uu .* S .- T(param_muuv) .* Vv

    return hcat(f1,f2,f3)
end

# ============================================================
# FULL OPERATOR
# ============================================================

function F_hat!(F_hat, U_hat, p)
    T = eltype(U_hat)

    # IMPORTANT: keep everything AD-consistent
    nu2_T = T.(nu2)
    lf = T.(laplace_factor)

    U_real = fourier_to_real(U_hat)
    Fnl_real = F_nonl_real(U_real)
    Fnl_hat = real_to_fourier(Fnl_real)

    for s in 1:N_species
        diffcoef = (s == bif_param) ? (p.diffcoef isa Number ? p.diffcoef : Float64(p.diffcoef)) : nu2_T[s]

        @views F_hat[:,s] .= lf .* U_hat[:,s] .* diffcoef .+ Fnl_hat[:,s]
    end

    return F_hat
end

function F_flat!(F_flat, U_flat, p)
    T = eltype(U_flat)

    Nk = size(C,2)
    U_hat = reshape(U_flat, Nk, N_species)

    F_hat = Array{T}(undef, Nk, N_species)

    F_hat!(F_hat, U_hat, p)

    F_flat .= vec(F_hat)
    return F_flat
end

# ============================================================
# INITIAL CONDITION
# ============================================================

U0 = [3.9354323319686637, 3.9354323319686637, 60.95051055800917]
U0_real = hcat(U0[1] .* ones(Nx),
               U0[2] .* ones(Nx),
               U0[3] .* ones(Nx))

U0_hat = real_to_fourier(U0_real)
sol0 = vec(U0_hat)

# ============================================================
# PARAMETER (CRITICAL FIX: MUST BE PURE Float64)
# ============================================================

par = (diffcoef = Float64(nu2[bif_param]),)

lens = (@optic _.diffcoef)

prob = BifurcationProblem(
    F_flat!,
    sol0,
    par,
    lens;
    record_from_solution = (x,p; k...) -> (
        nrm = norm(x),
        nrmReal = norm(x)
    ),
    plot_solution = (x,p; k...) -> begin
        T = eltype(x)
        C_T = T.(C)

        Nk = size(C,2)
        comp = x[(bif_param-1)*Nk+1 : bif_param*Nk]

        plot!(C_T * comp; k...)
    end
)

# ============================================================
# CONTINUATION SETTINGS
# ============================================================

opts = ContinuationPar(
    ds=1e-5,
    dsmax=5e-3,
    dsmin=1e-7,
    p_min=0.01,
    p_max=0.1,
    nev=2*N_fourier,
    n_inversion=6,
    detect_bifurcation=3,
    newton_options = NewtonPar(tol=1e-10, max_iterations=20),
    max_steps=200
)
# ============================================================
# RUN
# ============================================================

diagram = bifurcationdiagram(
    prob,
    PALC(),
    3,
    opts;
    bothside=true,
    plot=true
)

plot(diagram)