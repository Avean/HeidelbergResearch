using GLMakie
using DifferentialEquations
using SparseArrays
using LinearAlgebra
using Random

GLMakie.activate!()

# ============================================================
# 1D Laplacian with homogeneous Neumann boundary conditions
# ============================================================
function neumann_laplacian_1d(N, dx)
    L = spzeros(Float64, N, N)

    # left boundary: u_xx ≈ 2(u₂-u₁)/dx²
    L[1, 1] = -2.0
    L[1, 2] =  2.0

    # interior
    for i in 2:N-1
        L[i, i-1] =  1.0
        L[i, i]   = -2.0
        L[i, i+1] =  1.0
    end

    # right boundary: u_xx ≈ 2(u_{N-1}-u_N)/dx²
    L[N, N-1] =  2.0
    L[N, N]   = -2.0

    return L / dx^2
end


function fisher_kpp_demo(;
    N = 1000,
    D = 0.01,
    reltol = 1e-5,
    abstol = 1e-7,
    steps_per_frame = 5,
    dt_after_kick = 1e-8
)
    Random.seed!(1)

    # --------------------------------------------------------
    # Space grid
    # --------------------------------------------------------
    x = collect(range(0.0, 1.0; length=N))
    dx = x[2] - x[1]

    Lap = neumann_laplacian_1d(N, dx)
    u0 = (@. 0.2 + 0.05 * cos(2π * x)) + 0.02 * randn(N)

    # --------------------------------------------------------
    # RHS
    # --------------------------------------------------------
    function rhs!(du, u, p, t)
        mul!(du, Lap, u)
        @. du = D * du + u * (1.0 - u)
        return nothing
    end

    prob = ODEProblem(rhs!, u0, (0.0, Inf))

    # --------------------------------------------------------
    # Time-step controls & Time-shifting state
    # --------------------------------------------------------
    dt_choices = [1e-5, 1.0, 1e5, 1e-2]
    dtmax0 = 1e-2

    integrator_ref = Ref(init(
        prob,
        TRBDF2();
        adaptive = true,
        dt = min(1e-4, dtmax0),
        dtmax = dtmax0,
        reltol = reltol,
        abstol = abstol,
        save_everystep = false
    ))

    # Akumulator czasu – chroni przed utratą precyzji Float64
    time_offset = Ref(0.0)

    # --------------------------------------------------------
    # Makie observables
    # --------------------------------------------------------
    uobs = Observable(copy(integrator_ref[].u))
    tobs = Observable(integrator_ref[].t + time_offset[]) # UI widzi czas skumulowany
    dtmax_obs = Observable(dtmax0)
    dt_obs = Observable(integrator_ref[].dt)
    running_obs = Observable(false)
    step_counter = Observable(0) 

    title_obs = lift(tobs, dtmax_obs, dt_obs, running_obs, step_counter) do t, dtmax, dt, running, steps
        "Fisher-KPP | TRBDF2 | running = $(running) | t = $(round(t; digits=2)) | max dt = $(dtmax) | current dt = $(dt) | Steps: $(steps)"
    end

    fig = Figure(size = (1200, 760))
    ax = Axis(fig[1:14, 1], xlabel = "x", ylabel = "u", title = title_obs)
    lines!(ax, x, uobs; linewidth = 2)

    Label(fig[1, 2], "Wybierz max dt", tellwidth = false)

    # --------------------------------------------------------
    # Pomocnicze mechanizmy zarządzania czasem
    # --------------------------------------------------------
    simlock = ReentrantLock()

    function current_internal_dt()
        try return integrator_ref[].dt catch; return NaN end
    end

    # KLUCZOWA FUNKCJA: Przesunięcie czasu do zera bez straty stanu układu
    function shift_time_to_zero_if_needed!()
        # Jeśli wewnętrzny czas przekroczy 10,000, "zerujemy" zegar zachowując sumę w offset
        if integrator_ref[].t > 10000.0
            time_offset[] += integrator_ref[].t
            current_u = copy(integrator_ref[].u)
            current_dtmax = integrator_ref[].opts.dtmax
            
            new_prob = ODEProblem(rhs!, current_u, (0.0, Inf))
            integrator_ref[] = init(
                new_prob, TRBDF2(); 
                adaptive = true, dt = min(current_internal_dt(), current_dtmax), 
                dtmax = current_dtmax, reltol = reltol, abstol = abstol, save_everystep = false
            )
        end
    end

    function rescale_yaxis!()
        u = integrator_ref[].u
        umin, umax = minimum(u), maximum(u)
        margin = max(0.2 * (umax - umin), 0.1)
        ylims!(ax, umin - margin, umax + margin)
    end

    function refresh_plot!()
        uobs[] = copy(integrator_ref[].u)
        tobs[] = integrator_ref[].t + time_offset[] # Przekazujemy sumaryczny czas do wykresu
        dt_obs[] = current_internal_dt()
        dtmax_obs[] = integrator_ref[].opts.dtmax
        rescale_yaxis!()
    end

    function set_dtmax!(new_dtmax)
        if !(isfinite(new_dtmax) && new_dtmax > 0) return end
        lock(simlock)
        try
            integrator_ref[].opts.dtmax = new_dtmax
            cdmax = current_internal_dt()
            if isfinite(cdmax) && cdmax > new_dtmax
                set_proposed_dt!(integrator_ref[], new_dtmax)
            end
            refresh_plot!()
        finally unlock(simlock) end
    end

    # --------------------------------------------------------
    # RESTART PO PERTURBACJI (Z resetem lokalnego czasu)
    # --------------------------------------------------------
    function restart_after_manual_change!(unew)
        # Dodajemy dotychczasowy czas do offsetu i restartujemy solver od 0.0
        time_offset[] += integrator_ref[].t
        current_dtmax = max(1e-8, integrator_ref[].opts.dtmax)

        new_prob = ODEProblem(rhs!, unew, (0.0, Inf))

        integrator_ref[] = init(
            new_prob, TRBDF2();
            adaptive = true, dt = dt_after_kick, dtmax = current_dtmax,
            reltol = reltol, abstol = abstol, save_everystep = false
        )

        step_counter[] = 0
        refresh_plot!()
    end

    # --------------------------------------------------------
    # Kicks
    # --------------------------------------------------------
    function kick!(idx = :all; amount = 1.0)
        lock(simlock)
        try
            unew = copy(integrator_ref[].u)
            if idx === :all unew .+= amount else unew[idx] += amount end
            restart_after_manual_change!(unew)
        finally unlock(simlock) end
    end

    function sin_kick!(; amount = 1.0, mode = 1)
        lock(simlock)
        try
            unew = copy(integrator_ref[].u)
            @. unew += amount * sin(2π * mode * x)
            restart_after_manual_change!(unew)
        finally unlock(simlock) end
    end

    function sin_kick_clipped!(; amount = 1.0, mode = 1)
        lock(simlock)
        try
            unew = copy(integrator_ref[].u)
            @. unew += amount * sin(2π * mode * x) + 1.5
            @. unew = max(unew, 0.0)
            restart_after_manual_change!(unew)
        finally unlock(simlock) end
    end

    function step_once!()
        lock(simlock)
        try
            step!(integrator_ref[])
            step_counter[] += 1
            shift_time_to_zero_if_needed!()
            refresh_plot!()
        finally unlock(simlock) end
    end

    # --------------------------------------------------------
    # Pętla asynchroniczna z automatycznym shiftem czasu
    # --------------------------------------------------------
    task_ref = Ref{Union{Nothing, Task}}(nothing)

    function start!()
        if task_ref[] !== nothing && !istaskdone(task_ref[])
            running_obs[] = true
            return
        end
        running_obs[] = true

        task_ref[] = @async begin
            while running_obs[]
                lock(simlock)
                try
                    for _ in 1:steps_per_frame
                        step!(integrator_ref[])
                        step_counter[] += 1 
                    end
                    # Bezpiecznik sprawdzany co klatkę animacji
                    shift_time_to_zero_if_needed!()
                    refresh_plot!()
                catch err
                    @error "Błąd krytyczny w pętli integratora:" exception=(err, catch_backtrace())
                    running_obs[] = false
                finally unlock(simlock) end

                yield()
                sleep(0.001)
            end
        end
    end

    function stop!() running_obs[] = false end

    # --------------------------------------------------------
    # UI
    # --------------------------------------------------------
    for (i, dtchoice) in enumerate(dt_choices)
        b = Button(fig[i+1, 2], label = "max dt = $(dtchoice)", tellwidth = false)
        on(b.clicks) do _ set_dtmax!(dtchoice) end
    end

    Label(fig[6, 2], "Custom max dt:", tellwidth = false)
    dtbox = Textbox(fig[7, 2], placeholder = "np. 2e-3", stored_string = string(dtmax0), tellwidth = false)
    on(dtbox.stored_string) do s
        val = tryparse(Float64, s)
        if val !== nothing set_dtmax!(val) end
    end

    bstart, bstop, bone = Button(fig[9, 2], label="Start"), Button(fig[10, 2], label="Stop"), Button(fig[11, 2], label="Jeni krok")
    bkick = Button(fig[12, 2], label="kick: u += 1")
    bsinkick = Button(fig[13, 2], label="sin kick")
    bsinkickclip = Button(fig[14, 2], label="sin kick clipped")

    on(bstart.clicks)       do _ start!() end
    on(bstop.clicks)        do _ stop!() end
    on(bone.clicks)         do _ step_once!() end
    on(bkick.clicks)        do _ kick!() end
    on(bsinkick.clicks)     do _ sin_kick!() end
    on(bsinkickclip.clicks) do _ sin_kick_clipped!() end

    refresh_plot!()
    display(fig)
    start!()

    return (; fig, integrator_ref, step_counter)
end

ctrl = fisher_kpp_demo();