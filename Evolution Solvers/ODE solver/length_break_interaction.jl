# using Pkg
# Pkg.add("DifferentialEquations")
# Pkg.add("Plots")

using DifferentialEquations
using Roots
using Plots

# --- Parameters ---
const L = 1.0
const l_crit = 5.0 # Need l_crit^(-1) - C_out < 0
const b_crit = 1.0
const c_influx = 0.05
const M_start = 0.7
const C_out = 0.5

# --- Define l_des(t) implicitly ---
function compute_l_des(t)
    # f(l_des) = ((l_des - L) / l_des) * (1 / l_crit) +
    #            ((b_crit * c_influx * t + M_start) / l_des^2 + C_out)
    l_des_p = L / (2 - 2*l_crit*C_out) + sqrt(L^2 / (4*(1-l_crit*C_out)^2) - (b_crit + c_influx*t + M_start) / (l_crit^(-1) - C_out))
    # Choose a reasonable interval for root finding
    # result = find_zero(f, (0.01, 100.0), Bisection())  # safe range
    return l_des_p
end

# --- ODE System ---
function nonlinear_system!(du, u, p, t)
    b, l, M = u
    ldes = compute_l_des(t)
    du[1] = (- b + M*b)*(b_crit-b) #(ldes - l)
    term1 = (l - L) / (l * l_crit)
    term2 = M / l^2  # M = ((b_crit - b) * c_influx * t + M_start)
    du[2] = (term1 + term2 - C_out) * (l_crit - l)
    du[3] = c_influx*(b_crit-b) - M*b^3
end

# --- Initial Conditions and Time Span ---
u0 = [0.1, 1.2, M_start]  # initial values for b and l
tspan = (0.0, 100.0)  # from t = 0 to t = 10

# --- Define and solve the problem ---
prob = ODEProblem(nonlinear_system!, u0, tspan)
alg = Tsit5()
sol = solve(prob, alg)

# --- Plot ---
plot(sol, xlabel="Time", ylabel="Variables", label=["b(t)" "l(t)" "M(t)"], lw=2)
# ylims!(0,2)
