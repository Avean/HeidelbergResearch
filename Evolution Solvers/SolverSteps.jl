

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
    global DiffMat = undef
    global BC = undef
    global Order = undef
    global Type = undef

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

    function Choice(Scheme::String, BC::String, Order::String, NonFun::String)
        
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
        return SchemeF, f
    end

    function UpdateDiffMatrix()
        Solvers.DiffMat = CreateDiffMatrix(Sets.Par.Diff, DictLaplace[Solvers.Type], DictSinCosFun[Solvers.BC], DictSinCosEig[Solvers.BC], SimParam.dt);
        println("Diffusion matrix updated")
    end

    function run_simulation!(U0::VariablesVector,
                         Scheme::String,
                         BC::String,
                         Order::String,
                         Nonlinearity::String)

        # Globalisation in module
        
        Solvers.BC = BC;
        Solvers.Order = Order;
        Solvers.Type = BC*" "*Order;

        # --- Initialisation ---
        Fields    = fieldnames(VariablesVector)
        NFields   = 1:length(Fields)

        U1a = U0
        U1b = U0
        t   = 0.0
        V1  = [0.0]

        t1 = time()
        # --- Preparing time step functions ---
        SchemeF, FNonlinear = Choice(Scheme, BC, Order, Nonlinearity);
        
        
        UpdateDiffMatrix()

        println( time() - t1)
        
        # --- Main Simulation loop ---
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
                        Solvers.DiffMat[h],
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

