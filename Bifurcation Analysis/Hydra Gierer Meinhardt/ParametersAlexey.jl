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

includet("StructuresModule.jl")


using .Structures: Variables, Parameters

using .Nonlinearities

ν = [0.0, 3.8154e-05, 0.4433, 6.0713e-08, 0.0004];
β = [1.0629, 540.4003, 1.1596, 11.5964, 11.5964, 4.8254];

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


X0 = ones(5)
Xs = nlsolve((F,x)-> Fun!(F,x,Par), X0) # Trzeba to dodać bo Parametry jeszcze, funkcja anonimowa
Structures.StructPrint(Structures.Vec2Str(Xs.zero))