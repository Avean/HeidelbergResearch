module Viewer

using GLMakie
using Observables
using Printf
using Distributed
using ..SimParam          # dla L, N
using ..Struktury         # dla VariablesVector
using ..SharedState       # dla latest_state

export run_viewer

    function run_viewer()
        # 1. startowe wartości zanim solver zacznie słać prawdziwe dane
        U0  = VariablesVector(zeros(SimParam.N), zeros(SimParam.N))  # <- dostosuj do swoich pól
        V10 = [0.0]
        t0  = 0.0

        UObs  = Observable(U0)
        V1Obs = Observable(V10)
        tt    = Observable(t0)

        fig = Figure(title = "Hydra", resolution = (1920, 1080))
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
                                    "t = $(Printf.@sprintf("%0.1f", t_now))"
                end
            end
        end

        # --- dolne wykresy: variance i l-L ---
        tmin  = lift(t -> max(0, t - 30), tt)
        t1    = lift((t, tm) -> range(tm, t, 300), tt, tmin)
        tv    = lift((t, vhist) -> range(0, t, length(vhist)), tt, V1Obs)
        tgrow = lift(t -> (t ./ 10 .- floor.(t ./ 10)), t1)

        nvars = length(fieldnames(VariablesVector))

        axv = Axis(fig[nvars + 1, 1], xlabel = "Time", ylabel = "Variance")
        lines!(axv, tv, V1Obs)

        axt = Axis(fig[nvars + 2, 1],
                title = "Length l-L over time",
                xlabel = "Time",
                ylabel = "l-L")
        lines!(axt, t1, tgrow)

        on(tt) do ta
            xlims!(axt, max(0, ta - 30), ta)
            ylims!(axt, -0.1, 1.1)

            xlims!(axv, max(0, ta - 30), ta)
            vals = V1Obs[]
            ymin =minimum(vals) - 0.1*(maximum(vals) - minimum(vals)) 
            ymax = maximum(vals) + 0.1*(maximum(vals) - minimum(vals))
            ylims!(axv, ymin, ymax)
        end

        # 2. pętla odbiorcza 60 FPS
        @async begin
            fps = 60
            dt  = 1/fps
            while true
                # jeśli solver coś wysłał, zabierz najnowsze
                while isready(SharedState.latest_state)
                    U_now, V1_now, t_now = take!(SharedState.latest_state)
                    # ważne: kopiujemy tutaj, żeby Makie nie dostało obiektu,
                    # który solver jeszcze zmienia
                    UObs[]  = deepcopy(U_now)
                    V1Obs[] = deepcopy(V1_now)
                    tt[]    = t_now
                end
                sleep(dt)
            end
        end

        return nothing
    end
end