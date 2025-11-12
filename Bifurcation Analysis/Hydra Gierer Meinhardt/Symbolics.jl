
module SymbolicJacobian
    
    using ModelingToolkit


    export JacNonlinearity, FunNonlinearity

    @variables WntLoc DkkA WntDiff DkkC SD
    @parameters β1 β2 β3 β4 β5 β6

    F1 = β6 *SD / (1+DkkA) / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
    F2 = β1 / (1+ β4 * WntLoc) - DkkA
    F3 = β2 * WntLoc * SD - WntDiff
    F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
    F5 = WntLoc - SD


    F = [F1 F2 F3 F4 F5]

    
    Var = [WntLoc DkkA WntDiff DkkC SD];
    β = [β1 β2 β3 β4 β5 β6];

    FunNonlinearity, = ModelingToolkit.build_function(F,Var, β, expression = Val(false))
    
    H =  ModelingToolkit.jacobian(F,Var);
    JacNonlinearity,  = ModelingToolkit.build_function(H, Var, β,expression = Val(false));

end

