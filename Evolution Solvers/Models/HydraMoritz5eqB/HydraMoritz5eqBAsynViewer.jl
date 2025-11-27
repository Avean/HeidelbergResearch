module Viewer

using GLMakie
using Observables
using Printf
using Sockets
using Serialization
using Statistics
using DataStructures


# using Distributed
using ..SimParam          # dla L, N
using ..Struktury         # dla VariablesVector
using ..SharedState       # dla latest_state
using ..Nonlinearity
using ..WindowCreator
using ..Sets

# export run_viewer
export setup_viewer
export viewer_loop
export stop_simulation!



    fps = 60

    function ResetVariables!(U0,tdom, tval, V10, dt::Float64)


        tval[1:end] .= 0.0
        V10[1:end] .= 0.0
        tdom[1:end] .= 0.0
        U0  = Sets.Ini

        for i in 1:capacity(tdom)
            tdom[i] = (i-capacity(tdom))*dt
        end
    end

    function stop_simulation!(XVars)
        SharedState.pause_simulation[] = false

        SharedState.stop_simulation[] =true
        
        sleep(0.5)
        
        V10, tdom, tval, UObs, V1Obs, tt, ttv = XVars
        Viewer.ResetVariables!(UObs[],tt[], ttv[],V1Obs[], SimParam.dt)
        notify(V1Obs)
        notify(tt)
        notify(ttv)
    end


    function setup_viewer()

        U0  = Sets.Ini  
        THistory = 50; 

        tdom = CircularBuffer{Float64}(Int(floor(THistory*fps)))
        tval = CircularBuffer{Float64}(Int(floor(THistory*fps)))
        V10 = CircularBuffer{Float64}(Int(floor(THistory*fps)))

        fill!(tdom, 0.0)
        fill!(tval, 0.0)
        fill!(V10, 0.0)

        ResetVariables!(U0,tdom, tval, V10, SimParam.dt)

        UObs  = Observable(U0)
        V1Obs = Observable(V10)
        tt    = Observable(tdom)
        ttv   = Observable(tval)

        fig = Figure(title = "Hydra", resolution = (1020, 780))
        WindowCreator.show_fig!(WindowCreator.screen1, fig)

        Ω = range(0, SimParam.L, SimParam.N)
        
        ax = Vector(undef,fieldcount(Diffusions)) 


        # --- Upper plots for Variables ---
        for (i, fieldname) in enumerate(fieldnames(VariablesVector))
            Xi = lift(u -> getfield(u, fieldname), UObs)
            Xo = getfield(Sets.Ini, fieldname)

            ax[i] = Axis(fig[i, 1],
                    title = string(fieldname),
                    xlabel = "Ω",
                    ylabel = "Concentration")

            lines!(ax[i], Ω, Xi)
            lines!(ax[i], Ω, getfield(Sets.Ini,i), color = :red, linestyle = :dash)

            let ax_local = ax[i],
                X_local  = Xi,
                name_local = string(fieldname)

                on(tt) do t_now
                    data_now = X_local[]
                    ymin = min(minimum(data_now),minimum(Xo))
                    ymax = max(maximum(data_now),maximum(Xo))
                    Dist = ymax - ymin
                    # ymax = max(maximum(data_now) + 0.1, 2)
                    ylims!(ax_local, ymin-1e-8-0.1 * Dist, ymax+1e-8+0.1 * Dist)
                    ax_local.title[] = name_local * "\n" *
                                    "t = $(Printf.@sprintf("%0.1f", maximum(t_now)))"
                end
            end
        end

        XVars = (V10, tdom, tval, UObs, V1Obs, tt, ttv)
        return XVars
    end

    
    

    function viewer_loop!(X)
        V10, tdom, tval, UObs, V1Obs, tt, ttv = X
        
        time_start = time()
        display("Viewer started...")
        SharedState.stop_viewer[] = false
        while !SharedState.stop_viewer[]
            # println("Viewer loop time: ", time() - time_start)
            # time_start = time()
            # New frame:
            SharedState.request_frame[] = true

            # poczekaj aż symulacja ją wypełni
            # (bardzo krótko; symulacja biega w kółko i zauważy to prawie od razu)
            # proste aktywne czekanie:
            while SharedState.request_frame[]
                # jeszcze nie gotowe, daj innym wątkom czas
                yield()
            end

            # teraz frame_buffer[] powinno być ustawione
            time_start = time()
            snap = SharedState.frame_buffer[]
            if snap !== nothing
                U_now, t_now = snap

                # push!(V10, var(U_now.u))
                push!(tdom, t_now)
                # push!(tval, TimeSlope(Par,t_now))

                UObs[]  = U_now

                notify(V1Obs)
                notify(tt)
                notify(ttv)
            end

            sleep(1/fps)
        end
    end

    function server_loop!(Par::Parameters, X)
        V10, tdom, tval, UObs, V1Obs, tt, ttv = X
        fps = 60
        dt = 1/fps
        display("Viewer started...")
        while true
            # ta = time()
            snap = get_snapshot()
            if snap !== nothing
                U_now, t_now = snap

                push!(V10, var(U_now.u))
                push!(tdom, t_now)
                push!(tval, TimeSlope(Par,t_now))

                UObs[]  = U_now

                notify(V1Obs)
                notify(tt)
                notify(ttv)
            end
            # println("Viewer loop time: ", time() - ta)
            # sleep(dt)
        end

    end

    function get_snapshot(host::AbstractString = "127.0.0.1", port::Int = 20000)
        # println("Connecting to solver at $host:$port for snapshot...")
        sock = connect(host, port)
        # println("Connected, waiting for snapshot data...")
        snap = deserialize(sock)
        # println("Snapshot data received.")
        close(sock)
        return snap
    end
end