    using ModelingToolkit
    using Nemo
    using Symbolics

    
    function IsReal(x)
        if abs(imag(x)) < 1e-14 
            return real(x)
        else
            return X
        end
    end

    @variables WntLoc DkkA WntDiff DkkC SD
    
    # β1 = β[1]
    # β2 = β[2]
    # β3 = β[3]
    # β4 = β[4]
    # β5 = β[5]
    # β6 = β[6]


    β1 = 0.0;
    β2 = 1.5;
    β3 = 1.2;
    β4 = 0.0; 
    β5 = 1.5;
    β6 = 4.0;

    Var = [WntLoc, DkkA, WntDiff, DkkC, SD]

    # Nonliearity DkkA
    DkkA = β1 / (1+ β4 * WntLoc) *0
    DkkAV = ModelingToolkit.build_function(DkkA, Var,  expression=Val{false})
    
    
    # Nonliearity WntDiff
    WntDiff = β2 * WntLoc ^2
    WntDiffV = ModelingToolkit.build_function(WntDiff, Var,  expression=Val{false})
    
    
    # Nonliearity DkkC
    DkkC = WntDiff / (1 + β5 * WntLoc)
    DkkCV = ModelingToolkit.build_function(DkkC, Var,  expression=Val{false})
        
    # Nonliearity SD
    SD = WntDiff * 0.0
    SDV = ModelingToolkit.build_function(SD, Var,  expression=Val{false})

    



    
    F1 = β6 * WntDiff / (1+DkkC) / (1 + DkkA) / (1+ β3 * WntLoc) - WntLoc
    F1_frac = expand_derivatives(simplify(F1))
    num, den = numerator(F1_frac), denominator(F1_frac)

    # Coef = Symbolics.coefficients(num, WntLoc)
    X = (symbolic_solve(num~0, WntLoc))
    Y = Symbolics.symbolic_to_float.(X)
    Z = IsReal.(Y) |> collect |> x -> float.(sort(x, rev = true)) 
    
    U0 = []
    for i in Z
        push!(U0, [i; DkkAV(i); WntDiffV(i); DkkCV(i); SDV(i)])
    end

function pretty_H(H)
    L = latexify(H)  # albo bez env, jak wolisz
    println(L)                       # wypisz jako zwykły tekst
    return nothing                   # żeby REPL nie wyświetlał jeszcze raz stringa
end

H_sub = substitute(H, dict(WntDiff => WntDiff))