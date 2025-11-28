# α1 = b1 / c1
# α2 = b2 / c2
# α3 = b3 / c3
# α4 = b4 / c4
# α5 = b5 / c5

# β1 = α2*k1
# β2 = α3*α4*k2 / α5 
# β3 = k3 / α5
# β4 = k4 / α5
# β5 = k5 / α5
# β6 = α1 * α5

# ν2 = a2 / c2
# ν3 = a3 / c3
# ν4 = a4 / c4
# ν5 = a5 / c5


using Revise
using LinearAlgebra
using JLD2
using Base.Threads

includet("StructuresModule.jl")
includet("Symbolics.jl")


using .Structures: Variables, Parameters

using .Nonlinearities
using .SymbolicJacobian

ν = [0.0, 3.8154e-05, 0.4433, 6.0713e-08, 0.0004];
β = [1.0629, 540.4003, 1.1596, 11.5964, 11.5964, 4.8254];
# β2 = [1.0629, 540.4003, 1.1596, 11.5964, 11.5964, 4.8254];

# Var = Variables([1.0; 2.0],[1.0; 2.0],[1.0; 2.0],[1.0; 2.0],[1.0; 2.0])
# Var = Variables(1.0,1.0,1.0,1.0,1.0)

Par = Parameters(ν, β)

# X = F1(Var, Par)
# Structures.StructPrint(X)


###

using NLsolve


function Fun!(Val,x,Par)
    X = Structures.Vec2Str(x)
    Y = Nonlinearities.F1(X,Par)
    Val[:] = Structures.Str2Vec(Y)
    return nothing
end


# X0 = rand(5,1)
# Xs = nlsolve((F,x)-> Fun!(F,x,Par), X0) # Trzeba to dodać bo Parametry jeszcze, funkcja anonimowa
# Structures.StructPrint(Structures.Vec2Str(Xs.zero))
# Xs0 = Xs.zero




# function Fun2!(V,x,β)
#     V .= vec(FunNonlinearity(x, β))
#     return nothing
# end


X0 = rand(5,1)
XsSymbolic = nlsolve(x-> FunNonlinearity(x, β), 
                    x-> JacNonlinearity(x, β),
                    X0) 
Structures.StructPrint(Structures.Vec2Str(XsSymbolic.zero))

###

ZeroSol = [0.0; β[1]; 0.0; 0.0; 0.0];
@load "Points.jld2" A

    Eig1 = []
    for k in K
        EigVal = eigvals(JacNonlinearity(Var0,β) - Diagonal(ν) .*k^2.0.*π^2)
        Eig1 = [Eig1; maximum(real(EigVal))]
        # Eig0 = [Eig0; real.(EigVal[argmax((real.(EigVal)))])]
    end

# display(norm(FunNonlinearity(B, β)))

sim = Threads.@spawn for k in range(1,1e4)
    # for j in range(1,1e3)
        X0 = 10*rand(5)
        Xo = nlsolve((F,x)-> Fun2!(F,x,β), X0)
        
        
        if norm(Xo.zero - XsSymbolic.zero) >1e-6 && 
            norm(Xo.zero - ZeroSol) >1e-6 &&
            all([norm(Xo.zero - X) >1e-6 for X in A]) &&
            (norm(FunNonlinearity(Xo.zero, β))) < 1e-8
            break
        end
    if k % 1000 == 0
        println(k)
    end
end


# Structures.StructPrint(Structures.Vec2Str(Xo.zero))
# display(norm(FunNonlinearity(Xo.zero, β)))

# B = Xo.zero;
