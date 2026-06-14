# src/Spatial.jl

# ============================================================
# Spatial discretization utilities
# ============================================================


function make_grid_1d(N::Int; xmin::Float64 = 0.0, xmax::Float64 = 1.0)
    # Create a uniform one-dimensional spatial grid.
    #
    # Returns:
    #
    #     x  - vector of grid points
    #     dx - spatial mesh size

    N >= 2 || error("N must be at least 2.")

    x = collect(range(xmin, xmax; length = N))
    dx = x[2] - x[1]

    return x, dx
end


function neumann_laplacian_1d(N::Int, dx::Float64)
    # Build a sparse finite-difference matrix for the one-dimensional Laplacian
    # with homogeneous Neumann boundary conditions.
    #
    # Boundary condition:
    #
    #     u_x = 0
    #
    # on both endpoints.
    #
    # The boundary approximation is:
    #
    #     u_xx(x_1) ≈ 2(u_2 - u_1) / dx^2
    #     u_xx(x_N) ≈ 2(u_{N-1} - u_N) / dx^2

    N >= 2 || error("N must be at least 2.")
    dx > 0 || error("dx must be positive.")

    L = spzeros(Float64, N, N)

    # Left boundary
    L[1, 1] = -2.0
    L[1, 2] =  2.0

    # Interior points
    for i in 2:N-1
        L[i, i-1] =  1.0
        L[i, i]   = -2.0
        L[i, i+1] =  1.0
    end

    # Right boundary
    L[N, N-1] =  2.0
    L[N, N]   = -2.0

    return L / dx^2
end