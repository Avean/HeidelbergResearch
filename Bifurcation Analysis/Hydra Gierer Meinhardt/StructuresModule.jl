
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

module Nonlinearities
    using ..Structures: Variables, Parameters

    function F1(Var::Variables, Par::Parameters)
        VOut = Variables(
                        Par.β[6] .* Var.SD ./ (1 .+ Var.DkkA) ./ (1.0 .+ Var.DkkC) ./ (1.0 .+ Par.β[3]*Var.WntLOC) .- Var.WntLOC,
                        Par.β[1] ./ (1 .+ Par.β[4].*Var.WntLOC) .- Var.DkkA,
                        Par.β[2] .* Var.WntLOC.*Var.SD .- Var.WntDiff,
                        Var.WntDiff ./ ( 1.0 .+ Par.β[5] .* Var.WntLOC) .- Var.DkkC,
                        Var.WntLOC .- Var.SD,
                        ) 
        return VOut 
    end
end