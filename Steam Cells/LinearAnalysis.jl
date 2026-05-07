using ModelingToolkit
using Symbolics
using Symbolics: polynomialize
using SymbolicUtils
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)
using Plots


function GatherD(Dexp, A)
    Da = Differential(A)
    c0 = substitute(Dexp, A=>0)
    c1 = substitute(Da(Dexp) |> expand_derivatives, A=>0)
    c2 = substitute(Da(Da(Dexp))/2 |> expand_derivatives, A=>0)
    c3 = substitute(Da(Da(Da(Dexp)))/6 |> expand_derivatives, A=>0)

    return [c3, c2, c1, c0]    
end

@variables x y z 


f = [x^2, x+y]

J = ModelingToolkit.jacobian(f, [x, y])

J_fun = (build_function(J, x, y; expression=Val{false}))[1]


@variables A Q K H p0 r0 pa ra pq rq, β, b0
par = Dict(K=>1.01, H=>4.0, p0=>0.5, r0=>0.7 , pa=>0.1, ra=>0.2, pq=>0.3, rq=>0.4)

Up = pa * A + pq * Q
Vp = A + Q
p = p0 * (Up) / (H + Vp)

Ur = ra* A + rq* Q
Vr = A^2 + Q
r = r0 * (Vr) / (K + Vr)

B = b0 / (1.0 + β * Q)


L = A * Up * (K + Vr)
R = Q * Ur * (H + Vp)

Le = expand(L)
Re = expand(R)

Dexp = Le - Re

Coef = GatherD(Dexp, A)

C_num = substitute(Coef, par)

L_num = substitute(L, par) 
R_num = substitute(R, par) 

L_ex = expand(L_num)
R_ex = expand(R_num)

L_fun = (build_function(L_num, A, Q; expression=Val(false)))
R_fun = (build_function(R_num, A, Q; expression=Val(false)))

x = 0.0:0.1:10.0

y1 = L_fun.(x, 1.0)
y2 = R_fun.(x, 1.0)

plot(x, y1, label="L", xlabel="A", ylabel="Value")
plot!(x, y2, label="R")





