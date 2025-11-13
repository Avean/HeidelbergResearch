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

# export run_viewer
export setup_viewer
export viewer_loop
export stop_simulation!



    fps = 60

    function ResetVariables!(U0,tdom, tval, V10, dt::Float64)


        tval[1:end] .= 0.0
        V10[1:end] .= 0.0
        tdom[1:end] .= 0.0
        U0  = VariablesVector(zeros(SimParam.N), zeros(SimParam.N), zeros(SimParam.N))

        for i in 1:capacity(tdom)
            tdom[i] = (i-capacity(tdom))*dt
        end
    end

    function stop_simulation!(XVars, dt::Float64)
        SharedState.pause_simulation[] = false

        SharedState.stop_simulation[] =true

        V10, tdom, tval, UObs, V1Obs, tt, ttv = XVars
        Viewer.ResetVariables!(UObs[],tt[], ttv[],V1Obs[], dt)
        notify(V1Obs)
        notify(tt)
        notify(ttv)
    end


    function setup_viewer(Par::Parameters,dt::Float64)

        U0  = VariablesVector(zeros(SimParam.N), zeros(SimParam.N))  # <- dostosuj do swoich pól
        THistory = 50; 

        tdom = CircularBuffer{Float64}(Int(floor(THistory*fps)))
        tval = CircularBuffer{Float64}(Int(floor(THistory*fps)))
        V10 = CircularBuffer{Float64}(Int(floor(THistory*fps)))

        fill!(tdom, 0.0)
        fill!(tval, 0.0)
        fill!(V10, 0.0)

        ResetVariables!(U0,tdom, tval, V10, dt)

        UObs  = Observable(U0)
        V1Obs = Observable(V10)
        tt    = Observable(tdom)
        ttv   = Observable(tval)

        fig = Figure(title = "Hydra", resolution = (1020, 780))
        display(fig)

        Ω = range(0, SimParam.L, SimParam.N)
        Names = ["Morphogen Concentration u", "v"]

        # --- górne wykresy u, v ---
        for (i, fieldname) in enumerate(fieldnames(VariablesVector))
            Xi = lift(u -> getfield(u, fieldname), UObs)

            ax = Axis(fig[i, 1],
                    title = Names[i],
                    xlabel = "Ω",
                    ylabel = "Concentration")

            lines!(ax, Ω, Xi)

            let ax_local = ax,
                X_local  = Xi,
                name_local = Names[i]

                on(tt) do t_now
                    data_now = X_local[]
                    ymin = -0.1
                    ymax = max(maximum(data_now) + 0.1, 2)
                    ylims!(ax_local, ymin, ymax)
                    ax_local.title[] = name_local * "\n" *
                                    "t = $(Printf.@sprintf("%0.1f", maximum(t_now)))"
                end
            end
        end

        nvars = length(fieldnames(VariablesVector))

        axv = Axis(fig[nvars + 1, 1], xlabel = "Time", ylabel = "Variance")
        lines!(axv, tt, V1Obs)

        axt = Axis(fig[nvars + 2, 1],
                title = "Value of κ(l-L) over time",
                xlabel = "Time",
                ylabel = "κ(l-L)")
        lines!(axt, tt, ttv)
        ylims!(axt, -0.1, Par.Coef.lbreak + 0.1)

        on(tt) do tz
            ta = maximum(tz)
            xlims!(axt, max(0, ta - THistory), max(ta, 0.1))

            xlims!(axv, max(0, ta - THistory), max(ta,0.1))
            vals = V1Obs[]
            vals = vals[1:min(max(end-100,2),end)]
            ymin =minimum(vals) - 0.1*(maximum(vals) - minimum(vals))-1e-8 
            ymax = maximum(vals) + 0.1*(maximum(vals) - minimum(vals))+1e-8
            ylims!(axv, ymin, ymax)
        end

        XVars = (V10, tdom, tval, UObs, V1Obs, tt, ttv)
        return XVars
    end

    
    function setup_viewer2Graphs(Par::Parameters,dt::Float64)

        U0  = VariablesVector(zeros(SimParam.N), zeros(SimParam.N), zeros(SimParam.N))  # <- dostosuj do swoich pól
        THistory = 50; 

        tdom = CircularBuffer{Float64}(Int(floor(THistory*fps)))
        tval = CircularBuffer{Float64}(Int(floor(THistory*fps)))
        V10 = CircularBuffer{Float64}(Int(floor(THistory*fps)))

        fill!(tdom, 0.0)
        fill!(tval, 0.0)
        fill!(V10, 0.0)

        ResetVariables!(U0,tdom, tval, V10, dt)

        UObs  = Observable(U0)
        V1Obs = Observable(V10)
        tt    = Observable(tdom)
        ttv   = Observable(tval)

        fig = Figure(title = "Hydra", resolution = (1020, 780))
        display(fig)

        Ω = range(0, SimParam.L, SimParam.N)
        Names = ["u", "v", "w"]

        # --- górne wykresy u, v ---
        for (i, fieldname) in enumerate(fieldnames(VariablesVector))
        # i = 1
        # fieldname = :u    
        Xi = lift(u -> getfield(u, fieldname), UObs)

            ax = Axis(fig[i, 1],
                    title = Names[i],
                    xlabel = "Ω",
                    ylabel = "Concentration")

            lines!(ax, Ω, Xi)

            let ax_local = ax,
                X_local  = Xi,
                name_local = Names[i]

                on(tt) do t_now
                    data_now = X_local[]
                    ymin = -0.1
                    ymax = max(maximum(data_now) + 0.1, 2)
                    ylims!(ax_local, ymin, ymax)
                    ax_local.title[] = name_local * "\n" *
                                    "t = $(Printf.@sprintf("%0.1f", maximum(t_now)))"
                end
            end
        end

        nvars = 1

        # axv = Axis(fig[nvars + 1, 1], xlabel = "Time", ylabel = "Variance")
        # lines!(axv, tt, V1Obs)

        # axt = Axis(fig[nvars + 1, 1],
        #         title = "Value of κ(l-L) over time",
        #         xlabel = "Time",
        #         ylabel = "κ(l-L)")
        # lines!(axt, tt, ttv)
        # ylims!(axt, -0.1, Par.Coef.lbreak + 0.1)

        # on(tt) do tz
        #     ta = maximum(tz)
        #     xlims!(axt, max(0, ta - THistory), max(ta, 0.1))

        #     # xlims!(axv, max(0, ta - THistory), max(ta,0.1))
        #     vals = V1Obs[]
        #     vals = vals[1:min(max(end-100,2),end)]
        #     ymin =minimum(vals) - 0.1*(maximum(vals) - minimum(vals))-1e-8 
        #     ymax = maximum(vals) + 0.1*(maximum(vals) - minimum(vals))+1e-8
        #     # ylims!(axv, ymin, ymax)
        # end

        XVars = (V10, tdom, tval, UObs, V1Obs, tt, ttv)
        return XVars
    end

    
    function setup_viewer3Graphs(Par::Parameters,dt::Float64)

        U0  = VariablesVector(zeros(SimParam.N), zeros(SimParam.N))  # <- dostosuj do swoich pól
        THistory = 50; 

        tdom = CircularBuffer{Float64}(Int(floor(THistory*fps*2)))
        tval = CircularBuffer{Float64}(Int(floor(THistory*fps*2)))
        V10 = CircularBuffer{Float64}(Int(floor(THistory*fps*2)))

        fill!(tdom, 0.0)
        fill!(tval, 0.0)
        fill!(V10, 0.0)

        ResetVariables!(U0,tdom, tval, V10, dt)

        UObs  = Observable(U0)
        V1Obs = Observable(V10)
        tt    = Observable(tdom)
        ttv   = Observable(tval)

        fig = Figure(title = "Hydra", resolution = (1020, 780))
        display(fig)

        Ω = range(0, SimParam.L, SimParam.N)
        Names = ["Morphogen Concentration u", "v"]

        # --- górne wykresy u, v ---

        i = 1
        fieldname = fieldnames(VariablesVector)[i]
        Xi = lift(u -> getfield(u, fieldname), UObs)

            ax = Axis(fig[i, 1],
                    title = Names[i],
                    xlabel = "Ω",
                    ylabel = "Concentration")

            lines!(ax, Ω, Xi)

            let ax_local = ax,
                X_local  = Xi,
                name_local = Names[i]

                on(tt) do t_now
                    data_now = X_local[]
                    ymin = -0.1
                    ymax = max(maximum(data_now) + 0.1, 2)
                    ylims!(ax_local, ymin, ymax)
                    ax_local.title[] = name_local * "\n" *
                                    "t = $(Printf.@sprintf("%0.1f", maximum(t_now)))"
                end
            end
        

        nvars = 1

        axv = Axis(fig[nvars + 1, 1], xlabel = "Time", ylabel = "Variance")
        lines!(axv, tt, V1Obs)

        axt = Axis(fig[nvars + 2, 1],
                title = "Value of κ(l-L) over time",
                xlabel = "Time",
                ylabel = "κ(l-L)")
        lines!(axt, tt, ttv)

        on(tt) do tz
            ta = maximum(tz)
            xlims!(axt, max(0, ta - THistory), max(ta, 0.1))

            xlims!(axv, max(0, ta - THistory), max(ta,0.1))
            vals = V1Obs[]
            vals = vals[1:min(max(end-100,2),end)]
            ymin =minimum(vals) - 0.1*(maximum(vals) - minimum(vals))-1e-8 
            ymax = maximum(vals) + 0.1*(maximum(vals) - minimum(vals))+1e-8
            ylims!(axv, ymin, ymax)

            tvmax = maximum(ttv[])
            ylims!(axt, -0.1, max(Par.Coef.lbreak + 0.1,tvmax + 0.1)) 
        end

        XVars = (V10, tdom, tval, UObs, V1Obs, tt, ttv)
        return XVars
    end

    function viewer_loop!(Par::Parameters, X)
        V10, tdom, tval, UObs, V1Obs, tt, ttv = X
        
        time_start = time()
        display("Viewer started...")
        SharedState.stop_viewer[] = false



        while !SharedState.stop_viewer[]
            # println("Viewer loop time: ", time() - time_start)
            # time_start = time()
            # poproś o nową klatkę:
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

                push!(V10, var(U_now.u))
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