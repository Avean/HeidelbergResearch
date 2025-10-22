
module Solvers
    using ..SimParam
    using ..Struktury
    using ..LaplaceDiscretisation
    using ..Nonlinearity
    using ..Dictionaries
    using ..DiffMat
    
    using Plots
    using LinearSolve
    using LinearAlgebra
    using SparseArrays
    using Printf

    export Iteration

    # Solvers
    function ExpliciteEuler(U::Vector{Float64}, DiffMat::DiffusionMat, FU::Vector{Float64}, dt::Float64)
        return U + dt.*(DiffMat.LapMat*U ./ (SimParam.dx.^2) + FU);
    end

    function IMEX(U::Vector{Float64}, DiffMat::DiffusionMat, FU::Vector{Float64}, dt::Float64)
        return solve(LinearProblem(DiffMat.LapMatSparse,(U + dt .* FU)));
    end

    function SpectralSinCos(U::Vector{Float64}, DiffMat::DiffusionMat, FU::Vector{Float64}, dt::Float64)
        return DiffMat.Fun.Full'*((DiffMat.Fun.Trunc*(U + dt .* FU)).*DiffMat.Eig);
    end

    function Choice(Scheme::String, BC::String, Order::String, Par::Parameters, dt::Float64, NonFun::String)
        Type = BC*" "*Order;
        DiffMat = CreateDiffMatrix(Par.Diff, DictLaplace[Type], DictSinCosFun[BC], DictSinCosEig[BC], dt);


        if Scheme == "ExpliciteEuler"
            SchemeF = ExpliciteEuler;
        elseif Scheme == "IMEX"
            SchemeF = IMEX;
        elseif Scheme == "SpectralSinCos"
            SchemeF = SpectralSinCos;
        else
            error("Wrong Scheme");
        end
        f = NonlinearityFunction[NonFun];
        return SchemeF, DiffMat, f
    end

    #Choosing Solver and parameters
    
    function Iteration(U0::VariablesVector,  Par::Parameters, T::Float64, Scheme::String, BC::String, Order::String, dt::Float64, Nonlinearity::String)
        Fields = fieldnames(VariablesVector);
        FieldsNum = length(Fields);
        NFields = 1:FieldsNum;
        
        P = plot(layout = (FieldsNum,1));
        for i in NFields
            plot!(subplot = i, getfield(U0, Fields[i]), title = string(Fields[i]), ylims=(minimum(getfield(U0, Fields[i])) - 0.1, maximum(getfield(U0, Fields[i])) + 0.1));
        end
        display(P)

        sleep(1)
        U1a = U0;
        U1b = U0;


        
        SchemeF, DiffMat, FNonlinear = Choice(Scheme,BC,Order,Par,dt, Nonlinearity);
    
        t = 0.0;
        while t < T
            t1 = time()    
            for i = 1:100
                t = t + dt;
                U1b = FNonlinear(Par,U1a,t);
                U1a = VariablesVector(map(h -> 
                                        SchemeF(getfield(U1a,Fields[h]), 
                                                DiffMat[h],
                                                getfield(U1b,Fields[h]), 
                                                dt), 
                                        NFields)...);
            end

            display(time()-t1)
            display("t = $(Printf.@sprintf("%0.1e",t))")
            
            # P = plot(layout = (FieldsNum,1));
            P = plot(layout = (FieldsNum+1,1));
            for i in NFields
                X = getfield(U1a, Fields[i]);
                plot!(subplot = i, X, title = string(Fields[i]), label = "t = $(Printf.@sprintf("%0.1e",t))", ylims=(minimum([X;0.0]) - 0.1, maximum([X;2.0]) + 0.1));
            end

            tt = 0:0.001:t;
            plot!(subplot = 3, tt, (tt/10 - floor.(tt/10)), title = "Length l-L")
            display(P)

            sleep(0.01);
        end
        return U1a;
    end
end

module Extractor

    using ..Struktury

    export StructExtract

    function StructExtract(S::T) where {T}
        FNames = fieldnames(T);
        return map(h -> getfield(S, h), FNames);
    end

end