module LaplaceDiscretisation

    using  ..SimParam

    include("FillFunctions.jl")
    using .FillMatrix
    
    using SparseArrays
    using LinearAlgebra
    using Test

    export Lap, ISparse, Basis, SC, Eig#, LapSparse



    struct LaplacianMatrices
        Dir2::Matrix{Float64}
        Dir4::Matrix{Float64}
        Dir6::Matrix{Float64}
        Dir8::Matrix{Float64}
    
        Neu2::Matrix{Float64}
        Neu4::Matrix{Float64}
        Neu6::Matrix{Float64}
        Neu8::Matrix{Float64}
    
        Per2::Matrix{Float64}
        Per4::Matrix{Float64}
        Per6::Matrix{Float64}
        Per8::Matrix{Float64}
    end

    struct LaplacianMatricesSparse
        Dir2::SparseMatrixCSC{Float64, Int64}
        Dir4::SparseMatrixCSC{Float64, Int64}
        Dir6::SparseMatrixCSC{Float64, Int64}
        Dir8::SparseMatrixCSC{Float64, Int64}
    
        Neu2::SparseMatrixCSC{Float64, Int64}
        Neu4::SparseMatrixCSC{Float64, Int64}
        Neu6::SparseMatrixCSC{Float64, Int64}
        Neu8::SparseMatrixCSC{Float64, Int64}
    
        Per2::SparseMatrixCSC{Float64, Int64}
        Per4::SparseMatrixCSC{Float64, Int64}
        Per6::SparseMatrixCSC{Float64, Int64}
        Per8::SparseMatrixCSC{Float64, Int64}
    end

    struct Basis
        Full::Matrix{Float64}
        Trunc::Matrix{Float64}
    end

    struct SinCosBasis
        D::Basis; # Dirichlet
        N::Basis; # Neumann
        P::Basis; # Periodic
        H::Basis; # Half of periodic
    end

    struct  SinCosEig 
         S::Vector{Float64}; # Dirichlet
         C::Vector{Float64}; # Neumann
        SC::Vector{Float64}; # Periodic
         H::Vector{Float64}; # Half of periodic
    end    

    Dis2 = ([1.0, -2.0, 1.0]);
    Dis4 = ([-1/12, 4/3, -5/2, 4/3, -1/12]);
    Dis6 = ([1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90]);
    Dis8 = ([-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560]);

    
    Lap = LaplacianMatrices(FiniteDifferceDirichlet(Dis2),
                            FiniteDifferceDirichlet(Dis4),
                            FiniteDifferceDirichlet(Dis6),
                            FiniteDifferceDirichlet(Dis8),

                            FiniteDifferceNeumann(Dis2),
                            FiniteDifferceNeumann(Dis4),
                            FiniteDifferceNeumann(Dis6),
                            FiniteDifferceNeumann(Dis8),

                            FiniteDiffercePeriodic(Dis2),
                            FiniteDiffercePeriodic(Dis4),
                            FiniteDiffercePeriodic(Dis6),
                            FiniteDiffercePeriodic(Dis8));

    ISparse = sparse(I(SimParam.N));


    function SinCosBasisPeriodic()
        S = ones(2*SimParam.SicCosNodes + 1, SimParam.N)
        for i =1:SimParam.SicCosNodes
            S[2*i,:] = sqrt(2).*sin.(SimParam.x .*i ./SimParam.L .*2 .*pi);
            S[2*i+1,:] = sqrt(2).*cos.(SimParam.x.*i./SimParam.L .*2 .*pi);
        end
        return S
     end

     function SinCosBasisHalfPeriodic()
        S = ones(2*SimParam.SicCosNodes + 1, SimParam.N)
        for i = 1:2*SimParam.SicCosNodes
            S[i+1,:] = sqrt(2).*cos.(SimParam.x .*i ./SimParam.L .*2 .*pi);
        end
        return S
     end

     function SinCosBasisNeumann()
        S = ones(2*SimParam.SicCosNodes + 1, SimParam.N)
        for i = 1:2*SimParam.SicCosNodes
            S[i+1,:] = sqrt(2).*cos.(SimParam.x .*i ./SimParam.L .*pi);
        end
        return S
     end

     function SinCosBasisDirichlet()
        S = zeros(2*SimParam.SicCosNodes + 1, SimParam.N)
        for i =1:2*SimParam.SicCosNodes+1
            S[i,:] = sqrt(2).*sin.(SimParam.x .*i ./SimParam.L .*pi);
        end
        return S
     end

     function TruncatedBasis(S::Matrix{Float64})
        S = S / sqrt(SimParam.N - 1); #/ sqrt(SimParam.L);
        ST = deepcopy(S);
        ST[:,1] = ST[:,1]/2;
        ST[:,end] = ST[:,end]/2;
        return Basis(deepcopy(S), ST);
     end

    SC = SinCosBasis(TruncatedBasis(SinCosBasisDirichlet()), 
                     TruncatedBasis(SinCosBasisNeumann()), 
                     TruncatedBasis(SinCosBasisPeriodic()),
                     TruncatedBasis(SinCosBasisHalfPeriodic()));

    # Orthogonality tests. Checking if the basis is orthonormal 

    @test norm(SC.D.Trunc*SC.D.Full' - I(2*SimParam.SicCosNodes+1)) <= 1e-10
    @test norm(SC.N.Trunc*SC.N.Full' - I(2*SimParam.SicCosNodes+1)) <= 1e-10
    @test norm(SC.P.Trunc*SC.P.Full' - I(2*SimParam.SicCosNodes+1)) <= 1e-10

    
     SE = (1:1:(2*SimParam.SicCosNodes+1)).^2 .* pi^2 ./ SimParam.L^2;
     CE = [0; 1:1:(2*SimParam.SicCosNodes)].^2 .* pi^2 ./ SimParam.L^2;
    SCE = [0; kron(1:1:(SimParam.SicCosNodes),[1; 1])].^2 .* 4 .* pi^2 ./ SimParam.L^2;
    H = [0; 1:1:(2*SimParam.SicCosNodes)].^2 .* 4 .* pi^2 ./ SimParam.L^2;

    Eig = SinCosEig(SE,CE,SCE, H);
    
end

module Dictionaries
    using ..LaplaceDiscretisation

    export DictLaplace, DictLaplaceSparse, DictSinCosFun, DictSinCosEig

    DictLaplace = Dict(
                        "Dirichlet 2"  => Lap.Dir2,
                        "Dirichlet 4"  => Lap.Dir4,
                        "Dirichlet 6"  => Lap.Dir6,
                        "Dirichlet 8"  => Lap.Dir8,
                        "Neumann 2"    => Lap.Neu2,
                        "Neumann 4"    => Lap.Neu4,
                        "Neumann 6"    => Lap.Neu6,
                        "Neumann 8"    => Lap.Neu8,
                        "Periodic 2"   => Lap.Per2,
                        "Periodic 4"   => Lap.Per4,
                        "Periodic 6"   => Lap.Per6,
                        "Periodic 8"   => Lap.Per8
                        );

    # DictLaplaceSparse = Dict(
    #                         "Dirichlet 2"  => LapSparse.Dir2,
    #                         "Dirichlet 4"  => LapSparse.Dir4,
    #                         "Dirichlet 6"  => LapSparse.Dir6,
    #                         "Dirichlet 8"  => LapSparse.Dir8,
    #                         "Neumann 2"    => LapSparse.Neu2,
    #                         "Neumann 4"    => LapSparse.Neu4,
    #                         "Neumann 6"    => LapSparse.Neu6,
    #                         "Neumann 8"    => LapSparse.Neu8,
    #                         "Periodic 2"   => LapSparse.Per2,
    #                         "Periodic 4"   => LapSparse.Per4,
    #                         "Periodic 6"   => LapSparse.Per6,
    #                         "Periodic 8"   => LapSparse.Per8
    #                         );    

    DictSinCosFun = Dict(
                        "Dirichlet" => SC.D,
                        "Neumann"   => SC.N,
                        "Periodic"  => SC.P
                    );

    DictSinCosEig = Dict(
                        "Dirichlet" => Eig.S,
                        "Neumann"   => Eig.C,
                        "Periodic"  => Eig.SC
                    );
end

module DiffMat
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

    function FillMatrix(D::Float64,MatLap::Matrix{Float64}, FunVec::Basis, EigVec::Vector{Float64}, dt::Float64)
        
        return DiffusionMat(
            D .* MatLap, 
            ISparse - D .* dt ./ SimParam.dx ^2 .* sparse(MatLap),
            inv.(ones(SimParam.SicCosNodes*2 + 1) + D .* dt .* EigVec),
            FunVec);
    end

    function CreateDiffMatrix(DMat::Diffusions, LapMat::Matrix{Float64}, FunVec::Basis, EigVec::Vector{Float64}, dt::Float64)
        return (map(h-> FillMatrix(getfield(DMat,h),LapMat, FunVec, EigVec, dt),fieldnames(typeof(DMat))));
    end



end


