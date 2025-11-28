    
    using ModelingToolkit
    using NonlinearSolve


    @variables WntLoc DkkA WntDiff DkkC SD
    @parameters β1 β2 β3 β4 β5 β6

    F1 = β6 *SD / (1+DkkA) / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
    F2 = β1 / (1+ β4 * WntLoc) - DkkA
    F3 = β2 * WntLoc * SD - WntDiff
    # F3 = β2 * WntLoc + 0.0002 * β2 * SD - WntDiff
    F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
    F5 = WntLoc - SD


    F = [F1 F2 F3 F4 F5]

    
    Var = [WntLoc DkkA WntDiff DkkC SD];
    β = [β1 β2 β3 β4 β5 β6];

    FunNonlinearity, = ModelingToolkit.build_function(F,Var, β, expression = Val(false))
    
    H =  ModelingToolkit.jacobian(F,Var);
    JacNonlinearity,  = ModelingToolkit.build_function(H, Var, β,expression = Val(false));


    ν = [0.0, 3.8154e-05, 0.4433, 6.0713e-08, 0.0004];
    β = [1.0629, 540.4003, 1.1596, 11.5964, 11.5964, 4.8254];
    @load "SolP.jld2" Xs0


    @time Z = FunNonlinearity(U1, β)
    
    function Fun!(V,x,β)
        V .= vec(FunNonlinearity(x, β))
        return nothing
    end

    X0 = rand(5)
    Xs = nlsolve((x,F) -> Fun!(x,F,β), X0)
    JacNonlinearity(Xs.zero,β)


    function F(x::T) where T
        return x^2
    end






            function PrepareNonlinearity2(Fun, β, X0 = rand(fieldcount(VariablesVector)))
        F, J = Fun()
        U = SymbolicData.EvaluateConstant(Fun,β, X0)
        # if all(U .> 0.0 )
        #     return U, J(U,β)
        # end
    end