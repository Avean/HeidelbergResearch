# src/Spatial.jl

# ============================================================
# Spatial discretization utilities
# ============================================================

const SUPPORTED_BOUNDARY_CONDITIONS = (:neumann, :periodic)


function validate_boundary_condition(boundary_condition::Symbol)
    boundary_condition in SUPPORTED_BOUNDARY_CONDITIONS ||
        error("Unsupported boundary condition: $boundary_condition")

    return true
end


function boundary_condition_label(boundary_condition::Symbol)
    validate_boundary_condition(boundary_condition)

    if boundary_condition == :neumann
        return "Neumann"
    elseif boundary_condition == :periodic
        return "Periodic"
    end
end


function boundary_condition_from_label(label::String)
    if label == "Neumann"
        return :neumann
    elseif label == "Periodic"
        return :periodic
    else
        error("Unknown boundary condition label: $label")
    end
end


function boundary_condition_labels()
    return ["Neumann", "Periodic"]
end


function make_grid_1d(
    N::Int;
    xmin::Float64 = 0.0,
    xmax::Float64 = 1.0,
    boundary_condition::Symbol = :neumann,
)
    # Create a uniform one-dimensional spatial grid.
    #
    # For Neumann boundary conditions we include both endpoints:
    #
    #     x ∈ [xmin, xmax]
    #
    # For periodic boundary conditions we use a periodic grid without
    # duplicating the endpoint:
    #
    #     x ∈ [xmin, xmax)
    #
    # because xmin and xmax represent the same periodic point.

    validate_boundary_condition(boundary_condition)

    N >= 2 || error("N must be at least 2.")

    if boundary_condition == :neumann
        x = collect(range(xmin, xmax; length = N))
        dx = x[2] - x[1]

        return x, dx

    elseif boundary_condition == :periodic
        L = xmax - xmin
        L > 0 || error("xmax must be larger than xmin.")

        dx = L / N
        x = collect(xmin .+ dx .* (0:N-1))

        return x, dx
    end
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


function periodic_laplacian_1d(N::Int, dx::Float64)
    # Build a sparse finite-difference matrix for the one-dimensional Laplacian
    # with periodic boundary conditions.
    #
    # The grid is assumed to represent [xmin, xmax), so the last point is
    # connected back to the first point.

    N >= 3 || error("N must be at least 3 for periodic boundary conditions.")
    dx > 0 || error("dx must be positive.")

    L = spzeros(Float64, N, N)

    for i in 1:N
        im = i == 1 ? N : i - 1
        ip = i == N ? 1 : i + 1

        L[i, im] =  1.0
        L[i, i]  = -2.0
        L[i, ip] =  1.0
    end

    return L / dx^2
end


function laplacian_1d(
    N::Int,
    dx::Float64;
    boundary_condition::Symbol = :neumann,
)
    # Build a 1D Laplacian matrix for the selected boundary condition.

    validate_boundary_condition(boundary_condition)

    if boundary_condition == :neumann
        return neumann_laplacian_1d(N, dx)
    elseif boundary_condition == :periodic
        return periodic_laplacian_1d(N, dx)
    end
end