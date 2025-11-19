using Revise


includet("Symbolics.jl")
includet("ParametersAlexey.jl")

using LinearAlgebra
using Plots
using Statistics
using .SymbolicJacobian


Var0 = Xs.zero
β0 = β
ν0 = Diagonal(deepcopy(ν))


######## Eigenvalue closest to zero ###############

global A = []
K = 0:100
RangeΕ = [  
            range(10^(-0.64), 10^(-0.638), 300);
            range(10^(-1.227), 10^(-1.225), 300);
            range(10^(-1.554), 10^(-1.551), 300);
            range(10^(-1.62), 10^(-1.6), 300);
            range(10^(-1.768), 10^(-1.766), 300);
            range(10^(-1.92), 10^(-1.911), 1000);
            range(10^(-2.097), 10^(-2.0925), 300);
            range(10^(-2.0775), 10^(-2.075), 300);
            range(10^(-2.045), 10^(-2.04), 300);
            range(10^(-2.015), 10^(-2.013), 300);
            range(10^(-2.105), 10^(-2.1), 300);
        ]
# RangeΕ = range(0.007,0.44,1000)
# RangeΕ = [  range(10^(-1.92),10^(-1.9),1000);
#             range(10^(-2.13),1e-2,1000);
#             range(10^(-2),10^(-1.92),1000);
#             range(10^(-2.12),10^(-3),1000);
#             ]
# RangeΕ = range(10^(-2.13),1e-2,1000)



global Eig0 = []
global Colors = []

for ϵ in RangeΕ
    global Eig0 = []
    ν0[3,3] =  ϵ

    for k in K
        EigVal = eigvals(JacNonlinearity(Var0,β0) - ν0.*k^2.0.*π^2)
        Eig0 = [Eig0; real.(EigVal[argmin(abs.(real.(EigVal)))])]
    end
    # display(argmin(abs.(real.(Eig0))))
    push!(Colors, argmin(abs.(real.(Eig0))))
    # A = [A; Eig0[argmin(abs.(real.(Eig0)))]]
    push!(A,Eig0[argmin(abs.(real.(Eig0)))])
end
# scatter(K,real.(Eig0))

pCol = palette(:tab20)[1:14]
scatter(RangeΕ, real.(A), xscale = :log10, markersize = 1, markerstrokewidth = 0, marker_z = Colors, color = pCol,  colorbar_ticks = (1:14, string.(1:14)), colorbar = false)


p = plot()   # pusty wykres
yA = real.(A)
xA = RangeΕ

for k in 1:14
    idx = findall(Colors .== k)   # indeksy punktów o wartości z = k
    
    scatter!(p,
        xA[idx], yA[idx],
        markersize = 2,
        markerstrokewidth = 0,
        color = pCol[k],
        label = "k = $(k-1)",
        xscale = :log10

    )
end

ylims!(p,-0.0002,0.0002)
display(p)



mask = (-0.00001 .<yA .<0.00001)
yA1 = yA[mask]
xA1 = xA[mask]
Colors1 = Colors[mask]
scatter(Colors1,xA1,legend = false)

for i in unique(Colors1)
    println("Eigenode = ",i-1,"    Parameter = ", mean(xA1[Colors1.==i]),",")
end
1

##################################

#Max eigenvalue
# A = []
# for ϵ in RangeΕ
#     Eig0 = []
#     ν0[3,3] = ν[3]*ϵ

#     for k in R
#         EigVal = eigvals(JacNonlinearity(Var0,β0) + ν0.*k)
#         Eig0 = [Eig0; EigVal[argmax((real.(EigVal)))]]
#     end
#     A = [A; Eig0[argmax((real.(Eig0)))]]
# end
# scatter(K,real.(Eig0))
# scatter(RangeΕ, real.(A))

# ##
# k = 2
# Eig0 = []
# ν0[3,3] = 0.21

# EigVal = eigvals(JacNonlinearity(Var0,β0) - ν0.*k^2 * π^2)





# Eig0 = [Eig0; EigVal[argmax((real.(EigVal)))]]

# for k in R
#     EigVal = eigvals(JacNonlinearity(Var0,β0) + ν0.*k)
#     Eig0 = [Eig0; EigVal[argmax((real.(EigVal)))]]
# end

# scatter(K,real.(Eig0))