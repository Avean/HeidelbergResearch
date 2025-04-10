ModelName = "HydraDietmar";
# ModelName = "Test1";

######
###### Include the model files
######

includet(ModelName*"/"*ModelName*"Modules.jl")
include(ModelName*"/"*ModelName*"Variables.jl")

######
###### Solvers File
######

includet("DiffusionMatrices.jl")
includet("SolverSteps.jl")


module testAA
using ..Struktury
using ..Dictionaries
using ..SimParam
using ..LaplaceDiscretisation
using ..Sets
using ..Dictionaries

using SparseArrays

export DiffusionMat, CreateDiffMatrix 

struct DiffusionMat 
    LapMat::Matrix{Float64}
    LapMatSparse::SparseMatrixCSC{Float64, Int64}
    Eig::Vector{Float64}
    Fun::Basis
end

function CreateDiffStructure(input_struct::Diffusions)
    field_names = fieldnames(typeof(input_struct))
    
    new_struct_name = Symbol("DiffMats")
    

    struct_expr = Expr(:struct, false, new_struct_name, 
        Expr(:block, 
            map(field_names) do field_name
                Expr(:(::), field_name, :DiffusionMat)   
            end...
        )
    )
    eval(struct_expr)
            
    new_struct_type = eval(new_struct_name)
    return new_struct_type;

end

X = CreateDiffStructure(Sets.Set1.Diff)

function FillMatrix(D::Float64,MatLap::Matrix{Float64}, FunVec::Basis, EigVec::Vector{Float64}, dt::Float64)
    
    return DiffusionMat(
        D .* MatLap, 
        ISparse - D .* dt ./ SimParam.dx ^2 .* sparse(MatLap),
        inv.(ones(SimParam.SicCosNodes*2 + 1) + D .* dt .* EigVec),
        FunVec);
end



function CreateDiffMatrix(input_struct::Diffusions, LapMat::Matrix{Float64}, FunVec::Basis, EigVec::Vector{Float64}, dt::Float64)
    Smat = CreateDiffStructure(input_struct)
    display(fieldtypes((Smat)))

    field_names = fieldnames(typeof(input_struct))
    W =(map(h -> FillMatrix(getfield(input_struct,h), LapMat, FunVec, EigVec, dt),field_names));
    
    Smat(W...)
end
end

using ..testAA
using ..LaplaceDiscretisation

testAA.CreateDiffMatrix(Sets.Set1.Diff, Lap.Dir2, SC.N, Eig.S, 0.01)


X= plot(layout = (1,3));
plot!(subplot = 1, (1:10).^2);
plot!(subplot = 2, 1:10);

display(X)

Y= plot(layout = (1,2));

display(Y)

p1 = plot(rand(10), title="Pierwszy")
p2 = plot(rand(10), title="Drugi")
Plots.display(p1)
Plots.display(p2)



