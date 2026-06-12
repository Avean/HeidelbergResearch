module Viewer

using GLMakie
using Observables
using Printf
using Sockets
using Serialization
using Statistics
using DataStructures

using ..SimParam          # for L, N, dt
using ..Struktury         # for VariablesVector
using ..SharedState       # for shared simulation state
using ..Nonlinearity      # for TimeSlope with Sets.Par
using ..WindowCreator     # for show_fig!
using ..Sets              # for Sets.Ini, Sets.Par

export setup_viewer
export viewer_loop!
export stop_simulation!
export RecordAnimation
export server_loop!
export get_snapshot

# ============================================================================
# Constants and global state
# ============================================================================

const VIEWER_FPS = 60  # target frame rate for the viewer

# Global figure handle used e.g. by RecordAnimation
Fig = nothing


# ============================================================================
# Internal utilities
# ============================================================================

"""
    ResetVariables!(UObs, tdom, tval, V10, dt)

Reset time/history buffers to zero and reinitialize the time domain `tdom`
to a uniform grid going backwards in time with step `dt`. Also resets the
state observable `UObs` to `Sets.Ini`.

Parameters
----------
- `UObs` – observable holding the current state `U`
- `tdom` – circular buffer with time samples
- `tval` – circular buffer with κ(l-L) values (or similar)
- `V10`  – circular buffer with variance values
- `dt`   – time step used for constructing the time domain
"""
function ResetVariables!(
                            UObs::Observable,
                            tdom,
                            tval,
                            V10,
                            dt::Float64,
                        )
    tval[1:end] .= 0.0
    V10[1:end]  .= 0.0
    tdom[1:end] .= 0.0

    # Reset the state observable to the initial condition
    UObs[] = Sets.Ini
    notify(UObs)

    for i in 1:capacity(tdom)
        tdom[i] = (i - capacity(tdom)) * dt
    end
end


# ============================================================================
# Public API: stopping / setup
# ============================================================================

"""
    stop_simulation!(XVars)

Request stopping the simulation and the viewer, then reset viewer buffers.

`XVars` is the tuple:
    (V10, tdom, tval, UObs, V1Obs, tt, ttv)

This function:
- clears the `pause_simulation` flag,
- sets `stop_simulation` and `stop_viewer`,
- waits briefly to give other threads a chance to react,
- resets the viewer history buffers and the state observable,
- notifies the relevant observables.
"""
function stop_simulation!(XVars)
    
    SharedState.pause_simulation[] = false
    SharedState.stop_simulation[]  = true
    SharedState.stop_viewer[]      = true

    # Give other threads some time to react to stop signals.
    sleep(0.5)

    V10, tdom, tval, UObs, V1Obs, tt, ttv = XVars

    ResetVariables!(UObs, tt[], ttv[], V1Obs[], SimParam.dt)

    notify(V1Obs)
    notify(tt)
    notify(ttv)
end


"""
    setup_viewer() -> XVars

Initialize viewer buffers, observables and the Makie figure.

Returns
-------
A tuple
    (V10, tdom, tval, UObs, V1Obs, tt, ttv)
which should be passed to `viewer_loop!`, `server_loop!` and `stop_simulation!`.
"""
function setup_viewer()
    U0 = Sets.Ini
    THistory = 50.0

    nsteps = Int(floor(THistory * VIEWER_FPS))

    tdom = CircularBuffer{Float64}(nsteps)
    tval = CircularBuffer{Float64}(nsteps)
    V10  = CircularBuffer{Float64}(nsteps)

    fill!(tdom, 0.0)
    fill!(tval, 0.0)
    fill!(V10, 0.0)

    UObs  = Observable(U0)
    V1Obs = Observable(V10)
    tt    = Observable(tdom)
    ttv   = Observable(tval)

    # Reset buffers and state in a consistent way
    ResetVariables!(UObs, tdom, tval, V10, SimParam.dt)

    fig = Figure(title = "Hydra", resolution = (1020, 780))
    WindowCreator.show_fig!(WindowCreator.screen1, fig)

    Ω = range(0, SimParam.L, SimParam.N)

    # One axis for each variable in VariablesVector
    ax = Vector{Axis}(undef, fieldcount(VariablesVector))

    # --- Upper plots for variables ---
    for (i, fieldname) in enumerate(fieldnames(VariablesVector))
        Xi = lift(u -> getfield(u, fieldname), UObs)
        Xo = getfield(Sets.Ini, fieldname)

        ax[i] = Axis(
            fig[i, 1];
            title  = string(fieldname),
            xlabel = "Ω",
            ylabel = "Concentration",
        )

        lines!(ax[i], Ω, Xi)

        # Dynamically update axis limits and title based on the current data
        let ax_local   = ax[i],
            X_local    = Xi,
            name_local = string(fieldname),
            X_ref      = Xo

            on(tt) do t_now
                data_now = X_local[]
                ymin = min(minimum(data_now), minimum(X_ref))
                ymax = max(maximum(data_now), maximum(X_ref))
                dist = ymax - ymin

                ylow  = ymin - 1e-8 - 0.1 * dist
                yhigh = max(ymax + 1e-8 + 0.1 * dist, 0.0)

                ylims!(ax_local, ylow, yhigh)

                tmax = maximum(t_now)
                variance_val = var(data_now)

                ax_local.title[] = name_local * "\n" *
                    "t = $(Printf.@sprintf("%0.1f", tmax)),   " *
                    "Variance = $(Printf.@sprintf("%0.3g", variance_val))"
            end
        end
    end

    # Optionally enable extra plots (variance and κ(l-L) over time)
    # ExtraPlots(fig, V1Obs, tt, ttv, THistory)

    XVars = (V10, tdom, tval, UObs, V1Obs, tt, ttv)
    Viewer.Fig = fig

    return XVars
end


"""
    ExtraPlots(fig, V1Obs, tt, ttv, THistory)

Add additional time-series plots to the figure:

- variance of the first variable over time,
- κ(l-L) over time.

This function assumes:
- `tt`   is an observable for the time buffer,
- `V1Obs` is an observable for `V10` (variance history),
- `ttv`  is an observable for κ(l-L) values.
"""
function ExtraPlots(fig, V1Obs, tt, ttv, THistory)
    nvars = fieldcount(VariablesVector)

    axv = Axis(
        fig[nvars + 1, 1];
        xlabel = "Time",
        ylabel = "Variance",
    )
    lines!(axv, tt, V1Obs)

    axt = Axis(
        fig[nvars + 2, 1];
        title  = "Value of κ(l-L) over time",
        xlabel = "Time",
        ylabel = "κ(l-L)",
    )
    lines!(axt, tt, ttv)

    on(tt) do tz
        # Update x-limits to keep a sliding window of length THistory
        ta = maximum(tz)
        xlims!(axt, max(0, ta - THistory), max(ta, 0.1))
        xlims!(axv, max(0, ta - THistory), max(ta, 0.1))

        # Use a limited subset near the end to stabilize variance scaling
        vals = V1Obs[]
        vals = vals[1:min(max(end - 100, 2), end)]

        ymin = minimum(vals) - 0.1 * (maximum(vals) - minimum(vals)) - 1e-8
        ymax = maximum(vals) + 0.1 * (maximum(vals) - minimum(vals)) + 1e-8
        ylims!(axv, ymin, ymax)

        tvmax = maximum(ttv[])
        ylims!(axt, -0.1, max(Sets.Par.Coef.lbreak + 0.1, tvmax + 0.1))
    end
end


# ============================================================================
# Viewer loops
# ============================================================================

"""
    viewer_loop!(XVars)

Main viewer loop that runs in-process, communicating with the simulation
through `SharedState`.

`XVars` is the tuple:
    (V10, tdom, tval, UObs, V1Obs, tt, ttv)

Protocol
--------
On each iteration:
1. Sets `SharedState.request_frame[] = true` to ask the simulation for a new frame.
2. Actively waits (with `yield()`) until the simulation clears this flag and
   writes a snapshot into `SharedState.frame_buffer`.
3. If a snapshot is available, updates the viewer buffers and observables.
"""
function viewer_loop!(XVars)
    V10, tdom, tval, UObs, V1Obs, tt, ttv = XVars

    display("Viewer started...")
    SharedState.stop_viewer[] = false

    time_start = time()

    while !SharedState.stop_viewer[]
        time_start = time()

        # Request a new frame from the simulation
        SharedState.request_frame[] = true

        # Active wait: yield to give other threads time to respond
        while SharedState.request_frame[]
            yield()
        end

        # Now frame_buffer[] should be set by the simulation
        snap = SharedState.frame_buffer[]

        if snap !== nothing
            U_now, t_now = snap
            push!(V10, var(getfield(U_now, 1)))
            push!(tdom, t_now)
            # push!(tval, Nonlinearity.TimeSlope(Sets.Par, t_now))

            UObs[] = U_now

            notify(V1Obs)
            notify(tt)
            notify(ttv)
        end

        # Optional throttling:
        # sleep(1 / VIEWER_FPS)
    end
end


"""
    RecordAnimation(T, FileName, delay)

Record an animation of the current viewer figure over a time interval of length `T`
(seconds of viewer time), assuming that the observables are updated elsewhere
(e.g. by `viewer_loop!` or `server_loop!`).

Parameters
----------
- `T`        – length of the recording in seconds (in terms of frame count)
- `FileName` – output file name (e.g. `"movie.mp4"`)
- `delay`    – additional sleep time between frames (in real time), used to
               slow down recording if needed

Notes
-----
- Uses the global `Viewer.Fig` handle set in `setup_viewer`.
- Does not close the figure after recording (`close = false`).
"""
function RecordAnimation(T::Float64, FileName::AbstractString, delay::Float64)
    time_start = time()

    record(
        Viewer.Fig,
        FileName,
        1:Int(floor(T * VIEWER_FPS));
        close     = false,
        framerate = VIEWER_FPS,
    ) do i
        # We do not change any plot objects here; we simply wait for
        # the viewer loop to update the observables.
        display(time() - time_start)
        time_start = time()
        sleep(delay)
    end

    display("Recording saved to $FileName")
end


"""
    server_loop!(Par, XVars) - To be used later!

Viewer loop variant that receives snapshots from an external solver process
via TCP sockets instead of using `SharedState`.

Parameters
----------
- `Par`   – parameter object (type `Parameters` or similar), used by `TimeSlope`
- `XVars` – tuple (V10, tdom, tval, UObs, V1Obs, tt, ttv) created by `setup_viewer`
"""
function server_loop!(Par::Parameters, XVars)
    V10, tdom, tval, UObs, V1Obs, tt, ttv = XVars

    display("Viewer started...")

    while true
        snap = get_snapshot()
        if snap !== nothing
            U_now, t_now = snap

            push!(V10, var(U_now.u))
            push!(tdom, t_now)
            push!(tval, TimeSlope(Par, t_now))

            UObs[] = U_now

            notify(V1Obs)
            notify(tt)
            notify(ttv)
        end

        # Optional throttling:
        # fps = VIEWER_FPS
        # dt = 1 / fps
        # sleep(dt)
    end
end


# ============================================================================
# Socket communication
# ============================================================================

"""
    get_snapshot([host, port]) -> Any

Retrieve a snapshot from a running solver via TCP.

Parameters
----------
- `host` – hostname or IP address (default `"127.0.0.1"`)
- `port` – TCP port (default `20000`)

Returns
-------
Deserialized object sent by the solver, typically a tuple `(U_now, t_now)`.

Notes
-----
- The socket is closed in a `finally` block to ensure it is always released,
  even if deserialization fails.
"""
function get_snapshot(host::AbstractString = "127.0.0.1", port::Int = 20000)
    sock = connect(host, port)
    try
        snap = deserialize(sock)
        return snap
    finally
        close(sock)
    end
end

end
