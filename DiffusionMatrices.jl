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
        D::Basis;
        N::Basis;
        P::Basis;
    end

    struct  SinCosEig 
         S::Vector{Float64};  
         C::Vector{Float64};
        SC::Vector{Float64};
    end    

    Dis2 = ([1.0, -2.0, 1.0]);
    Dis4 = ([-1/12, 4/3, -5/2, 4/3, -1/12]);
    Dis6 = ([1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90]);
    Dis8 = ([-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560]);

    # function FiniteDiffercePeriodic(Var::Vector{Float64})
    #     L = zeros(SimParam.N, SimParam.N)
    #     k =  floor(Int, length(Var)/2 - 0.5);
    #     for i = 1:SimParam.N
    #         for j = -k:1:k
    #             L[i,mod(i-j-1,SimParam.N)+1] = Var[k+1-j];
    #         end
    #     end

    #     return L
    # end

    # function FiniteDifferceNeumann(Var::Vector{Float64})
    #     L = zeros(SimParam.N, SimParam.N)
    #     k =  floor(Int, length(Var)/2 - 0.5);
    #     for i = 1:SimParam.N
    #         L[i,i] = Var[k+1];
    #         for j = 1:k
    #             L[i,i-j>0 ? i-j : 1 ] = Var[k+1-j];
    #             L[i,i+j <SimParam.N ? i+j : SimParam.N] = Var[k+1+j];
    #         end
    #     end
    #     for i=1:k
    #         L[i, 1] = sum(Var[1:k+2-i])
    #         L[end-i+1, end] = sum(Var[k+i:end])
    #     end
    #     return L
    # end

    # function FiniteDifferceDirichlet(Var::Vector{Float64})
    #     L = zeros(SimParam.N, SimParam.N)
    #     k =  floor(Int, length(Var)/2 - 0.5);
    #     for i = 1:SimParam.N
    #         for j = -min(k,i-1):1:0
    #             L[i+j,i] = Var[k+1+j];
    #         end
    #         for j = 0:1:min(k,SimParam.N-i)
    #             L[i+j,i] = Var[k+1+j];
    #         end
    #     end
    #     return L
    # end
    
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

    # LapSparse = LaplacianMatricesSparse(sparse(FiniteDifferceDirichlet(Dis2)),
                                        # sparse(FiniteDifferceDirichlet(Dis4)),
                                        # sparse(FiniteDifferceDirichlet(Dis6)),
                                        # sparse(FiniteDifferceDirichlet(Dis8)),

                                        # sparse(FiniteDifferceNeumann(Dis2)),
                                        # sparse(FiniteDifferceNeumann(Dis4)),
                                        # sparse(FiniteDifferceNeumann(Dis6)),
                                        # sparse(FiniteDifferceNeumann(Dis8)),

                                        # sparse(FiniteDiffercePeriodic(Dis2)),
                                        # sparse(FiniteDiffercePeriodic(Dis4)),
                                        # sparse(FiniteDiffercePeriodic(Dis6)),
                                        # sparse(FiniteDiffercePeriodic(Dis8)));
                


    function SinCosBasisPeriodic()
        S = ones(2*SimParam.SicCosNodes + 1, SimParam.N)
        for i =1:SimParam.SicCosNodes
            S[2*i,:] = sqrt(2).*sin.(SimParam.x .*i ./SimParam.L .*2 .*pi);
            S[2*i+1,:] = sqrt(2).*cos.(SimParam.x.*i./SimParam.L .*2 .*pi);
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
                     TruncatedBasis(SinCosBasisPeriodic()));

    # Orthogonality tests. Checking if the basis is orthonormal 

    @test norm(SC.D.Trunc*SC.D.Full' - I(2*SimParam.SicCosNodes+1)) <= 1e-10
    @test norm(SC.N.Trunc*SC.N.Full' - I(2*SimParam.SicCosNodes+1)) <= 1e-10
    @test norm(SC.P.Trunc*SC.P.Full' - I(2*SimParam.SicCosNodes+1)) <= 1e-10

    
     SE = (1:1:(2*SimParam.SicCosNodes+1)).^2 .* pi^2 ./ SimParam.L^2;
     CE = [0; 1:1:(2*SimParam.SicCosNodes)].^2 .* pi^2 ./ SimParam.L^2;
    SCE = [0; kron(1:1:(SimParam.SicCosNodes),[1; 1])].^2 .* 4 .* pi^2 ./ SimParam.L^2;

    Eig = SinCosEig(SE,CE,SCE);
    
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


