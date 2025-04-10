using ..Struktury
using ..Dictionaries
using ..SimParam
using ..LaplaceDiscretisation
using ..Sets
using ..Dictionaries

using SparseArrays

export DiffusionMat, CreateDiffMatrix 
export CreateDiffStructure

function CreateDiffStructure(input_struct, type)
    field_names = fieldnames(input_struct)
    
    new_struct_name = Symbol("Diff"*string(type))
    
    struct_expr = Expr(:struct, false, new_struct_name, 
        Expr(:block, 
            map(field_names) do field_name
                Expr(:(::), field_name, :($(type)))   
            end...
        )
    )
    display(struct_expr)
    # Ewaluacja wyrażenia w celu stworzenia struktury
    eval(struct_expr)
    
    # Nowy typ struktury
   new_struct_type = eval(new_struct_name)
   return new_struct_type
end

T1 = CreateDiffStructure(Diffusions, Matrix{Float64});
T2 = CreateDiffStructure(Diffusions, SparseMatrixCSC{Float64});
T3 = CreateDiffStructure(Diffusions, Matrix{Float64});

display(fieldtypes(T1))
display(fieldtypes(T2))
display(fieldtypes(T3))   

struct DiffusionMat 
    LapMat::T1
    LapMatSparse::T2
    Eig::T3
end

X = DiffusionMat(T1([1;;]), 
                 T2(sparse([1;;])), 
                 T3([1;;]));

                #  MatLap = Lap.Dir2;
                #  EigVec = Eig.S;
                #  dt = 0.01;
# function CreateDiffMatrix(MatLap::Matrix{Float64}, EigVec::Vector{Float64}, dt::Float64)

    Fields = fieldnames(Diffusions)

    X1 = CreateDiffStructure(Diffusions, Matrix{Float64});
    DiffLap = X1(map(D -> 
                    D .* MatLap, 
                    Fields)...);

    X2 = CreateDiffStructure(Diffusions, SparseMatrixCSC{Float64});
    DiffLapSparse = X2(map(D -> 
                    ISparse - D .* dt ./ SimParam.dx ^2 .* sparse(MatLap),
                    Fields)...);

    X3 = CreateDiffStructure(Diffusions, Vector{Float64});    
    DiffEig = X3(map(D -> 
    inv.(ones(SimParam.SicCosNodes*2 + 1) + D .* dt .* EigVec), 
                    Fields)...);

    # return DiffusionMat(DiffLap, DiffLapSparse, DiffEig);

# end

CreateDiffMatrix(Lap.Dir2, Eig.S, 0.01);