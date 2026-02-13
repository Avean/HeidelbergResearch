    
    using ModelingToolkit
    using JLD2
    # using NonlinearSolve


    @variables WntLoc DkkA WntDiff DkkC SD
    @parameters β1 β2 β3 β4 β5 β6

    # F1 = β6 *SD / (1+DkkA) / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
    # F2 = β1 / (1+ β4 * WntLoc) - DkkA
    # F3 = β2 * WntLoc * SD - WntDiff
    # # F3 = β2 * WntLoc + 0.0002 * β2 * SD - WntDiff
    # F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
    # F5 = WntLoc - SD

            # ############# Nonlinearities ############
            # F1 = β6 * WntDiff / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
            # F2 = 0
            # F3 = β2 * WntLoc^2  - WntDiff
            # F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
            # F5 = 0           
            # ########################################


            # ############# Nonlinearities ############
            # # F1 = β6 * SD / (1+DkkC) / (1+ β3 * WntLoc) / (1+DkkA) - WntLoc
            # F1 = β6 * SD / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
            # # F2 = β1 / (1+ β4 * WntLoc) - DkkA
            # F2 = 0
            # F3 = β2 * WntLoc   - WntDiff
            # F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
            # F5 = WntLoc - SD           
            # ########################################

            ############# Nonlinearities ############
            F1 = β6 * WntDiff / (1+DkkC) / (1+ β3 * WntLoc) / (1+DkkA) - WntLoc
            # F1 = β6 * SD / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
            F2 = β1 / (1+ β4 * WntLoc) - DkkA
            # F2 = 0
            F3 = β2 * WntLoc * SD   - WntDiff
            F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
            F5 = WntLoc - SD           
            ########################################



                        ############# Nonlinearities ############
            # F1 = β6 * SD / (1+DkkC) / (1+ β3 * WntLoc) / (1+DkkA) - WntLoc
            # F1 = β6 * WntLoc / (1+DkkC) / (1+ β3 * WntLoc)  - WntLoc
            # # F2 = β1 / (1+ β4 * WntLoc) - DkkA
            # F2 = 0
            # F3 = β2 * WntLoc   - WntDiff
            # F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
            # F5 = 0           
            ########################################




    F = [F1 F2 F3 F4 F5]

    
    Var = [WntLoc DkkA WntDiff DkkC SD];
    β = [β1 β2 β3 β4 β5 β6];

    FunNonlinearity, = ModelingToolkit.build_function(F,Var, β, expression = Val(false))
    
    H =  ModelingToolkit.jacobian(F,Var);
    JacNonlinearity,  = ModelingToolkit.build_function(H, Var, β,expression = Val(false));


    ν = [0.0, 3.8154e-05, 0.4433, 6.0713e-08, 0.0004];
    β = [1.0629, 540.4003, 1.1596, 11.5964, 11.5964, 4.8254];
    @load "SolP.jld2" Xs0


    # @time Z = FunNonlinearity(U1, β)
    
    # function Fun!(V,x,β)
    #     V .= vec(FunNonlinearity(x, β))
    #     return nothing
    # end

    # X0 = rand(5)
    # Xs = nlsolve((x,F) -> Fun!(x,F,β), X0)
    # JacNonlinearity(Xs.zero,β)


    # function F(x::T) where T
    #     return x^2
    # end





using Latexify

"""
Zwraca ładny kod LaTeX dla zadanej macierzy H::Matrix{Num},
z podmianą nazw zmiennych:
  WntDiff -> W_D
  WntLoc  -> W_L
  DkkC/ DKKC -> C
  β2, β3, β5, β6 -> \beta_{2}, \beta_{3}, ...
i opakowaniem w \begin{align*} ... \end{align*}.
"""


function matrix_to_pretty_latex3(H)
    # 1. Generujemy surowy LaTeX
    lat = latexify(H)
    s = String(lat)

    # 2. USUWAMY \mathtt
    # Regex wyłapuje zawartość wewnątrz \mathtt{...} i zostawia tylko środek.
    # Wzorzec r"\\mathtt\{([^{}]*)\}" szuka \mathtt{cokolwiek_bez_nawiasow}
    s = replace(s, r"\\mathtt\{([^{}]*)\}" => s"\1")
    
    # Alternatywa (brutalna): po prostu usuń komendę, zostawiając nawiasy {} (LaTeX to zignoruje)
    # s = replace(s, "\\mathtt" => "") 

    # 3. Naprawa pmatrix (nawiasy okrągłe zamiast kwadratowych)
    s = replace(s, r"\\begin{equation\*?}" => "")
    s = replace(s, r"\\end{equation\*?}" => "")
    s = replace(s, "\\left[" => "")
    s = replace(s, "\\right]" => "")
    s = replace(s, r"\\begin{array}\{.*?\}" => "\\begin{pmatrix}")
    s = replace(s, "\\end{array}" => "\\end{pmatrix}")

    # 4. Podmiana nazw na symbole matematyczne
    # Ponieważ usunęliśmy \mathtt, teraz mamy czyste nazwy (lub z nawiasami klamrowymi z latexify)
    replacements = [
        "WntDiff"   => "W_D",
        "WntLoc"    => "W_L",
        "DkkC"      => "C",
        "DKKC"      => "C",
        "SD"       => "S",
        "DkkA"      => "A",
        # Latexify czasem robi {\beta}6, czasem β6, więc warto obsłużyć oba przypadki:
        r"\{\\beta\}6" => "\\beta_{6}",
        r"\{\\beta\}5" => "\\beta_{5}",
        r"\{\\beta\}3" => "\\beta_{3}",
        r"\{\\beta\}2" => "\\beta_{2}",
        r"\{\\beta\}1" => "\\beta_{1}",
        r"\{\\beta\}4" => "\\beta_{4}",
        "β6"        => "\\beta_{6}",
        "β5"        => "\\beta_{5}",
        "β3"        => "\\beta_{3}",
        "β2"        => "\\beta_{2}",
        "β1"        => "\\beta_{1}",
        "β4"        => "\\beta_{4}",
    ]

    for (old, new) in replacements
        s = replace(s, old => new)
    end

    return "\\begin{align*}\n" * strip(s) * "\n\\end{align*}"
end

println(matrix_to_pretty_latex3(H))