# using Pkg
# Pkg.add("DifferentialEquations")
# Pkg.add("Plots")

using DifferentialEquations
using Roots
using Plots

# --- Parameters ---
const L = 1.0
const l_crit = 5.0 # Need l_crit^(-1) - C_out < 0 (not anymore as we do not use l_des)
const b_crit = 1.0
const c_influx = 0.05
const M_start = 0.5
const C_out = 1.0

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
# function nonlinear_system!(du, u, p, t)
#     b, l, M = u
#     # ldes = compute_l_des(t)
#     du[1] = (- b + M*b)*(b_crit-b) #(ldes - l)
#     term1 = (l - L) / (l * l_crit)
#     term2 = M / l^2  # M = ((b_crit - b) * c_influx * t + M_start)
#     du[2] = (-term1 + term2 - C_out) * (l_crit - l) #-term1 or +term1?
#     du[3] = c_influx*(b_crit-b)*l - M*b^3 #produces oscillations: c_influx*M - M*b^3
#     # Need *l in first term of du[3] in order to obtain different speed of oscillations
#     # when C_out changed. Compare Fig.2 E in Ferenc 2021.
# end

function nonlinear_system!(du, u, p, t)
    b, l, M, conc = u
    # ldes = compute_l_des(t)
    conc_in = M /l^2
    du[1] = (- b + conc_in*b)*(b_crit-b) #(ldes - l)
    press = (l - L) / (l * l_crit)
    du[2] = (-press + conc_in - C_out) * (l_crit - l)
    # du[3] = (c_influx*conc_in - 2*conc_in^2*du[2]*l - 2*conc_in*b^3*l^3) # equation for conc_in
    du[3] = c_influx - b^3*conc_in  # equation for M
    du[4] = (du[3]*conc - 2*conc^2*du[2]*l) / M # equation for conc insides, we have conc = conc_in = M/l^2.
end

# --- Initial Conditions and Time Span ---
u0 = [0.01, L, M_start, M_start/L.^2]  # initial values for b and l and M (and conc_in: M_start/l^2)
tspan = (0.0, 200.0)  # from t = 0 to t = 10

# --- Define and solve the problem ---
prob = ODEProblem(nonlinear_system!, u0, tspan)
# alg = Tsit5()
sol = solve(prob)

# --- Plot ---
plot(sol, xlabel="Time", ylabel="Variables", label=["b(t)" "l(t)" "M(t)" "conc(t)"], lw=2)
l = sol[2,:]
M = sol[3,:]
concentration = M./l.^2
plot!(sol.t, concentration, label="cont_in(t)")
# ylims!(0,2)
