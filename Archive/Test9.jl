using BifurcationKit, Plots

function Fsl(X, p)
    (;r, μ, ν, c3) = p
    u, v = X
    ua = u^2 + v^2
    [
        r * u - ν * v - ua * (c3 * u - μ * v)
        r * v + ν * u - ua * (c3 * v + μ * u)
    ]
end



X = [0.1, 0.1];
p = (r = 1.0, μ = 0.1, ν = 0.1, c3 = 0.1);
Fsl(X, p)


par_sl = (r = 0.1, μ = 0., ν = 1.0, c3 = 1.0)
u0 = zeros(2)
prob = BifurcationProblem(Fsl, u0, par_sl, (@optic _.r))
br = continuation(prob, PALC(), ContinuationPar(), bothside = true)

br_po = continuation(br, 2, ContinuationPar(),
        PeriodicOrbitOCollProblem(20, 5)
        )
        br_po[1]

        plot(br, br_po, branchlabel = ["equilibria", "periodic orbits"])