using BifurcationKit, Plots
F(x, p) = @. (p[1] + x - x^3/3)*x
prob = BifurcationProblem(F, [-1.0], [-1.0], 1;
    record_from_solution = (x,p; k...) -> x[1])


opts_br = ContinuationPar(
    p_max = 1.0,
    p_min = -1.0,
    n_inversion = 4)

    diagram = bifurcationdiagram(prob, PALC(),
	# very important parameter. This specifies the maximum amount of recursion
	# when computing the bifurcation diagram. It means we allow computing branches of branches
	# at most in the present case.
	2,
	opts_br,
	)

# br = continuation(prob, PALC(tangent = Bordered()), ContinuationPar(p_min = -1., p_max = 1.))
plot(diagram)

