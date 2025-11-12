
includet("Symbolics.jl")

using LinearAlgebra
using Plots
using .SymbolicJacobian


Var0 = Xs.zero
β0 = β
ν0 = Diagonal(deepcopy(ν))


######## Eigenvalue closest to zero ###############

A = []
K = 0:50
R = [-π^2 * k^2 for k in K]
RangeΕ = range(-0.4,4,100)

for ϵ in RangeΕ
    Eig0 = []
    ν0[3,3] = ν[3] +  ϵ

    for k in R
        EigVal = eigvals(JacNonlinearity(VarO,βO) + ν0.*k)
        Eig0 = [Eig0; real.(EigVal[argmin(abs.(real.(EigVal)))])]
    end
    A = [A; Eig0[argmin(abs.(real.(Eig0)))]]
end
# scatter(K,real.(Eig0))
scatter(RangeΕ, real.(A))


##################################

#Max eigenvalue
A = []
for ϵ in RangeΕ
    Eig0 = []
    ν0[3,3] = ν[3] +  ϵ

    for k in R
        EigVal = eigvals(JacNonlinearity(VarO,βO) + ν0.*k)
        Eig0 = [Eig0; EigVal[argmax((real.(EigVal)))]]
    end
    A = [A; Eig0[argmax((real.(Eig0)))]]
end
# scatter(K,real.(Eig0))
scatter(RangeΕ, real.(A))

##
k = 0
Eig0 = []

EigVal = eigvals(JacNonlinearity(VarO,βO) - ν0.*k^2 * π^2)
Eig0 = [Eig0; EigVal[argmax((real.(EigVal)))]]

for k in R
    EigVal = eigvals(JacNonlinearity(VarO,βO) + ν0.*k)
    Eig0 = [Eig0; EigVal[argmax((real.(EigVal)))]]
end

scatter(K,real.(Eig0))