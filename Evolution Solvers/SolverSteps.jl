

module SharedState
    export request_frame, frame_buffer

    const request_frame = Ref(false)
    const frame_buffer  = Ref{Any}(nothing)
    const stop_simulation = Ref(false)
    const pause_simulation = Ref(false)
    const stop_viewer = Ref(false)

end

module Solvers
    using ..SimParam
    using ..Struktury
    using ..LaplaceDiscretisation
    using ..Nonlinearity
    using ..Dictionaries
    using ..DiffMat
    using ..SharedState
    using ..Sets
    
    # using Plots

    # using GLMakie
    # using Observables

    using LinearSolve
    using LinearAlgebra
    using SparseArrays
    using Printf
    using Statistics

    using Sockets
    using Serialization

    export Iteration

    const last_snapshot = Ref{Any}(nothing)

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
    
    # function Iteration(U0::VariablesVector,  Par::Parameters, T::Float64, Scheme::String, BC::String, Order::String, dt::Float64, Nonlinearity::String)
        
    #     Fields = fieldnames(VariablesVector)
    #     FieldsNum = length(Fields);
    #     NFields = 1:FieldsNum

    #     sleep(1)
    #     U1a = U0;
    #     U1b = U0;

    #     t = 0.0;
    #     V1 = [0.0]

        
    #     # GL Makie        
    #     # UObs, V1Obs, tt = DynamicPlotGLMakie(U1a, V1, t)

    #     SchemeF, DiffMat, FNonlinear = Choice(Scheme,BC,Order,Par,dt, Nonlinearity);
    
    #     while t < T
    #         t1 = time()    
    #             t = t + dt;
    #             U1b = FNonlinear(Par,U1a,t);
    #             U1a = VariablesVector(map(h -> 
    #                                     SchemeF(getfield(U1a,Fields[h]), 
    #                                             DiffMat[h],
    #                                             getfield(U1b,Fields[h]), 
    #                                             dt), 
    #                                     NFields)...);

    #         push!(V1, var(U1a.u))

    #                 # --- obsługa prośby o klatkę ---
    #         if SharedState.request_frame[]
    #             # kopiujemy snapshot żeby viewer dostał stabilne dane
    #             SharedState.frame_buffer[] = (deepcopy(U1a), deepcopy(V1), t)

    #             # sygnał "gotowe"
    #             SharedState.request_frame[] = false
    #         end

    #         # GL Makie
    #         # UObs[] = U1a
    #         # V1Obs[]=V1
    #         # tt[] = t

            
    #         # display(typeof(SharedState.latest_state))
    #         # spróbuj wyczyścić kanał, jeśli coś tam jeszcze siedzi
    #         # if isready(SharedState.latest_state)
    #             # take!(SharedState.latest_state)
    #         # end
    #         # teraz włóż najnowszy stan
    #         # put!(SharedState.latest_state, (U1a, V1, t))
            
    #         # V1Obs=push!(V1, var(U1a.u))
    #         # DynamicPlotGR(U1a, t, FieldsNum, Fields, V1Obs);
            
    #         display(time()-t1)
    #         # display("t = $(Printf.@sprintf("%0.1e",t))")    
            
    #     end

    #     return U1a;
    # end

    function run_simulation!(U0::VariablesVector,
                         Scheme::String,
                         BC::String,
                         Order::String,
                         Nonlinearity::String)

        # --- inicjalizacja pól ---
        Fields    = fieldnames(VariablesVector)
        NFields   = 1:length(Fields)

        U1a = Sets.Ini
        U1b = Sets.Ini
        t   = 0.0
        V1  = [0.0]

        # --- przygotowanie funkcji kroku czasowego ---
        SchemeF, DiffMat, FNonlinear =
            Choice(Scheme, BC, Order, Sets.Par, SimParam.dt, Nonlinearity)

        # --- główna pętla symulacji ---
        SharedState.stop_simulation[] = false
        while !SharedState.stop_simulation[]
            
            if SharedState.pause_simulation[]
                yield()
                continue
            end
            t1 = time()    
            # krok czasowy
            t += SimParam.dt
            U1b = FNonlinear(Sets.Par, U1a, t)
            U1a = VariablesVector(map(h ->
                SchemeF(getfield(U1a, Fields[h]),
                        DiffMat[h],
                        getfield(U1b, Fields[h]),
                        SimParam.dt),
                NFields)...)

            # zapisz np. wariancję jako prosty przykład obserwacji
            # push!(V1, var(U1a.u))

            # --- sprawdź, czy viewer prosi o nową klatkę ---
            if SharedState.request_frame[]
                # kopiujemy stan, żeby viewer dostał spójne dane
                SharedState.frame_buffer[] = (deepcopy(U1a), t)

                # flaga zresetowana -> viewer wie, że dane gotowe
                SharedState.request_frame[] = false
            end

            # display(time()-t1)
            last_snapshot[] = (deepcopy(U1a), t)
            # sleep(5)
            # oddaj schedulerowi czas na inne wątki (np. viewer)
            yield()
        end
    end



    function snapshot_server!(port::Int = 20000)
        println("Starting snapshot server on port $port")
        server = listen(port)
        println("Snapshot server listening on port $port")

        @async begin
            while true
                sock = accept(server)
                @async begin
                    try
                        # println("Client connected for snapshot")
                        snap = last_snapshot[]
                        # println("snap = ", snap)

                        if snap === nothing
                            println("Warning: last_snapshot is nothing, sending placeholder")
                            serialize(sock, (:no_data,))
                        else    
                            # println("Sending snapshot data")
                            serialize(sock, snap)
                            # println("Snapshot sent")
                        end
                    catch err
                        @error "handler error" err
                    finally
                        close(sock)
                    end
                end
            end
        end

        return nothing
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

