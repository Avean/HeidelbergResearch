using SparseArrays

N = 100

U = rand(N)

function F(U)
    return U.*(1.0.-U)
end

D = spdiagm(0 => -2*ones(N), 1 => ones(N-1), -1 => ones(N-1))
D[1,1] = -1
D[end,end] = -1


FU = spdiagm(0 => F(U))

