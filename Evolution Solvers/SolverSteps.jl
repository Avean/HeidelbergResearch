module SharedState
    using Distributed
    export latest_state
    const latest_state = RemoteChannel(() -> Channel{Any}(1))
end

module Solvers
    using ..SimParam
    using ..Struktury
    using ..LaplaceDiscretisation
    using ..Nonlinearity
    using ..Dictionaries
    using ..DiffMat
    using ..SharedState
    
    # using Plots
    using GLMakie
    using Observables

    using LinearSolve
    using LinearAlgebra
    using SparseArrays
    using Printf
    using Statistics

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

    function DynamicPlotGR(U1a::VariablesVector, t::Float64, FieldsNum::Int64, Fields, V1::Vector{Float64})
            P = plot(layout = (FieldsNum+2,1));
            for i in 1:FieldsNum
                X = getfield(U1a, Fields[i])
                ymin = min(minimum(X), 0.0) - 0.1
                ymax = max(maximum(X), 2.0) + 0.1
                plot!( X;
                    subplot = i,
                    title   = string(Fields[i]),
                    label   = "t = $(Printf.@sprintf("%0.1e", t))",
                    ylims   = (ymin, ymax))
            end

            tt = (max(0,t-30)):0.1:t;
            push!(V1, var(U1a.u))

            plot!( tt, (tt ./ 10 .- floor.(tt ./ 10)),subplot = 3, title = "Length l-L", label=false)
            plot!( range(0,t,length(V1)), V1;subplot = 4, title = "Variance", xlabel="iteracja", label=false, xlim = (max(0,t-30), t))
            display(P)     
    end


    function DynamicPlotGLMakie(U1a::VariablesVector, V1::Vector{Float64}, t::Float64)
        UObs = Observable(U1a)
        V1Obs = Observable(V1)
        tt = Observable(t)
        
        X = []

        fig = Figure(title = "Hydra", resolution = (1920, 1080))
        display(fig)


        Ω = range(0, SimParam.L, SimParam.N)
        Names = ["Morphogen Concentration u", "v"]
        
        for (i, j) in enumerate(fieldnames(VariablesVector))
            X = lift(x-> getfield(x, j), UObs)
            ax = Axis(fig[i, 1], title = Names[i], xlabel = "Ω", ylabel = "Concentration")
            
            lines!(ax, Ω, X)

            let X_local = X, ax_local = ax
            on(tt) do Y 
                Names = ["Morphogen Concentration u", "v"]
                data_now = X_local[]
                ymin =  - 0.1
                ymax = max(maximum(data_now) + 0.1,2)               
               ylims!(ax_local, ymin, ymax)
               ax_local.title[] =Names[i] * "\n" *"t = $(Printf.@sprintf("%0.1f", Y))"
            end
            end
        end

        tmin = lift(t -> max(0, t - 30), tt)
        t0 = lift((t,V1Obs) -> range(0,max(t,1),max(length(V1),2)), tt,V1Obs)

        
        t1 = lift((t,tm) ->range(tm,t,300),tt,tmin)
        tv = lift((t,V1Obs) -> range(0,t,length(V1)), tt,V1Obs)
        tgrow = lift(t -> (t ./ 10 .- floor.(t ./ 10)), t1)
        
        axv = Axis(fig[length(fieldnames(VariablesVector))+1, 1], xlabel = "Time", ylabel = "Variance")
        lines!(axv,tv, V1Obs)

        axt = Axis(fig[length(fieldnames(VariablesVector))+2, 1], title = "Length l-L over time", xlabel = "Time", ylabel = "l-L")
        lines!(axt, t1, tgrow)
        
        on(tt) do ta
            xlims!(axt, max(0, ta - 30), ta)
            ylims!(axt, -0.1, 1.1)

            xlims!(axv, max(0, ta - 30), ta)
            # autolimits!(axv)
            ylims!(axv, -0.1*maximum(V1Obs[])-0.1*(maximum(V1Obs[]) - minimum(V1Obs[])), 1.1*maximum(V1Obs[]))
        end
        return UObs, V1Obs, tt
    end

    #Choosing Solver and parameters
    
    function Iteration(U0::VariablesVector,  Par::Parameters, T::Float64, Scheme::String, BC::String, Order::String, dt::Float64, Nonlinearity::String)
        
        Fields = fieldnames(VariablesVector)
        FieldsNum = length(Fields);
        NFields = 1:FieldsNum

        sleep(1)
        U1a = U0;
        U1b = U0;

        t = 0.0;
        V1 = [0.0]

        
        # GL Makie        
        # UObs, V1Obs, tt = DynamicPlotGLMakie(U1a, V1, t)

        SchemeF, DiffMat, FNonlinear = Choice(Scheme,BC,Order,Par,dt, Nonlinearity);
    
        while t < T
            t1 = time()    
                t = t + dt;
                U1b = FNonlinear(Par,U1a,t);
                U1a = VariablesVector(map(h -> 
                                        SchemeF(getfield(U1a,Fields[h]), 
                                                DiffMat[h],
                                                getfield(U1b,Fields[h]), 
                                                dt), 
                                        NFields)...);

            push!(V1, var(U1a.u))

            # GL Makie
            # UObs[] = U1a
            # V1Obs[]=V1
            # tt[] = t

            
            # display(typeof(SharedState.latest_state))
            # spróbuj wyczyścić kanał, jeśli coś tam jeszcze siedzi
            if isready(SharedState.latest_state)
                take!(SharedState.latest_state)
            end
            # teraz włóż najnowszy stan
            put!(SharedState.latest_state, (U1a, V1, t))
            
            # V1Obs=push!(V1, var(U1a.u))
            # DynamicPlotGR(U1a, t, FieldsNum, Fields, V1Obs);
            
            display(time()-t1)
            # display("t = $(Printf.@sprintf("%0.1e",t))")    
            
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

