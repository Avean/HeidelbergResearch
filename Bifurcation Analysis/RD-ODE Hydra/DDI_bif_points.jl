using LinearAlgebra
using NLsolve

# Define the nonlinearites
beta = [1.0629 540.4003 1.1596 11.5964 11.5964 4.8254];
F_nonl(Wl, A, Wd, C, S) = [beta[6]*S ./ ((1.0 .+ A).*(1.0 .+ C).*(1.0 .+ beta[3]*Wl)) - Wl,
                           beta[1] ./ (1.0 .+ beta[4]*Wl) - A,
                           beta[2]*Wl.*S - Wd,
                           Wd ./ (1.0 .+ beta[5]*Wl) - C,
                           Wl - S
                          ]

# Non-trivial constant steady state:
U0 = [0.08167547924306925, 0.54587711524379, 3.6049476660047937, 1.8514050545898464, 0.08167547924306925];

# Jacobian:
t2 = beta[3]*U0[1];
t3 = beta[5]*U0[1];
t4 = U0[2]+1.0;
t5 = U0[4]+1.0;
t6 = t2+1.0;
t7 = t3+1.0;
t8 = 1.0/t4;
t9 = 1.0/t5;
t10 = 1.0/t6;
jacSym = reshape([-beta[3]*beta[6]*U0[5]*t8*t9*t10^2-1.0, -beta[1]*beta[4]*1.0/(beta[4]*U0[1]+1.0)^2, beta[2]*U0[5],
                  -beta[5]*1.0/t7^2*U0[3], 1.0, -beta[6]*U0[5]*t8^2*t9*t10, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0,
                  1.0/t7, 0.0, -beta[6]*U0[5]*t8*t9^2*t10, 0.0, 0.0, -1.0, 0.0, beta[6]*t8*t9*t10, 0.0, beta[2]*U0[1], 0.0, -1.0],
                5,5);

eigvals(jacSym)

# Diffusion coefficients
d = 0.0088;
nu = [0.0, 3.8154e-05, d, 6.0713e-08, 0.0004];

k=8;
D = diagm(nu);
eigvals(jacSym - k^2*pi^2*D)


function disprel(wavenumb, diffcoef)
    det(jacSym - wavenumb^2*pi^2*diagm([0.0, 3.8154e-05, diffcoef, 6.0713e-08, 0.0004]))
end

function disprel_fix_k(diffcoef)
    disprel(8, diffcoef)    
end

nlsolve(disprel_fix_k, [0.9])

Jac = [[-1.08652  -0.0528344   0.0       -0.0286439   1.0],
  [-3.25103  -1.0         0.0        0.0         0.0],
  [44.1375    0.0        -1.0        0.0        44.1375],
 [-11.0262    0.0         0.513573  -1.0         0.0],
   [1.0       0.0         0.0        0.0        -1.0]]
