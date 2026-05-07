
# x = range(0,20,100);

# z = plot(x,Set.Coef.c .* x)
# plot!(x,Set.Coef.b.^2 .*(x .+ Set.Coef.a) .^2);
# display(z)
# V = 15.7;
# U = Set.Coef.b* (Set.Coef.a +V)
# U = sqrt(V.*Set.Coef.c)

using Symbolics, Nemo, Printf

@variables u, v, a, b, c

W = Symbolics.jacobian([u^2/(v+a)-b*u, u^2-c*v], [u, v])

# u = b*(a+v)

expr = b^2*(v+a)^2 - c*v
V0 = simplify(symbolic_solve(expr, v))
V0 = simplify(V0.*a./a)
U0 = simplify(b.*(V0.+a))

U0n = substitute(U0, Dict(a => Set.Coef.a, b => Set.Coef.b, c => Set.Coef.c))
V0n = substitute(V0, Dict(a => Set.Coef.a, b => Set.Coef.b, c => Set.Coef.c))

J1 = simplify(substitute(W, Dict(u => U0n[1], v => V0n[1], a => Set.Coef.a, b => Set.Coef.b, c => Set.Coef.c)))
J2 = simplify(substitute(W, Dict(u => U0n[2], v => V0n[2], a => Set.Coef.a, b => Set.Coef.b, c => Set.Coef.c)))

@printf("Punkt 1: (%1.3f, %1.3f). Jacobian 1 wynosi: \n", U0n[1], V0n[1])
display(J1)
@printf("Wyznacznik %0.3f:\n", det(J1))
@printf("\n")
@printf("Punkt 2: (%1.3f, %1.3f). Jacobian 1 wynosi:\n", U0n[2], V0n[2])
display(J2)
@printf("Wyznacznik %0.3f:\n", det(J2))

# D1 = 0.0001;
# D2 = 0.002;
k = 3;
JK = J1 - k.^2 .*pi.^2 .* [D1 0; 0 D2];
display(det(JK))
display(JK)

L = 10.0;
J = 0.0;
dL = 0.05;
WS = [];
for k=J:dL:L
    JK = J1 - k.^2 .*pi.^2 .*4 .* [D1 0; 0 D2];
    push!(WS,maximum(real(eigvals(JK))))
    # display(eigvals(JK))
end


display(plot((J:dL:L),WS))
plot!([J, L],[0.0, 0.0],linecolor = "black")