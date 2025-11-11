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


module Structures

    struct Parameters{T1,T2}
        ν::T1
        β::T2
    end

    struct Variables{T}
        WntLOC::T
        DkkA::T
        WntDiff::T
        DkkC::T
        SD::T
    end

    function StructPrint(obj)
        T = typeof(obj)
        for name in fieldnames(T)
            println("$(name): ", getfield(obj, name))
        end
    end

    function Str2Vec(obj::T) where T
        return [getfield(obj, names) for names in fieldnames(T) ]
    end

    function Vec2Str(obj::Vector{T}) where T
        return Variables([obj[i] for i in 1:length(fieldnames(Variables))]...)
    end

end

using .Structures: Variables, Parameters

nu = [0.0 3.8154e-05 0.4433 6.0713e-08 0.0004];
beta = [1.0629 540.4003 1.1596 11.5964 11.5964 4.8254];

Var = Variables(1.0,1.0,1.0,1.0,1.0)
# Var = Variables([1.0; 2.0],[1.0; 2.0],[1.0; 2.0],[1.0; 2.0],[1.0; 2.0])

Par = Parameters(nu, beta)

function Nonlinearity(Var::Variables, Par::Parameters)
    VOut = Variables(
                    Par.β[6] .* Var.SD ./ (1 .+ Var.DkkA) ./ (1.0 .+ Var.DkkC) ./ (1.0 .+ Par.β[3]*Var.WntLOC) .- Var.WntLOC,
                    Par.β[1] ./ (1 .+ Par.β[4].*Var.WntLOC) .- Var.DkkA,
                    Par.β[2] .* Var.WntLOC.*Var.SD .- Var.WntDiff,
                    Var.WntDiff ./ ( 1.0 .+ Par.β[5] .* Var.WntLOC) .- Var.DkkC,
                    Var.WntLOC .- Var.SD,
                    ) 
    return VOut 
end

# X = Nonlinearity(Var, Par)
# Structures.StructPrint(X)


###

using NLsolve


function F!(Val,x,Par)
    X = Structures.Vec2Str(x)
    Y = Nonlinearity(X,Par)
    Val[:] = Structures.Str2Vec(Y)
    return nothing
end


X0 = ones(5)
Xs = nlsolve((F,x)-> F!(F,x,Par), X0) # Trzeba to dodać bo Parametry jeszcze, funkcja anonimowa
