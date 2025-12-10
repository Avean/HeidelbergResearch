    
    using ModelingToolkit
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


            ############# Nonlinearities ############
            F1 = β6 * WntDiff / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
            F2 = 0
            F3 = β2 * WntLoc^2  - WntDiff
            F4 = WntDiff / (1 + β5 * WntLoc)^2 - DkkC
            F5 = 0           
            ########################################

            F = [F1 F2 F3 F4 F5]


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
function matrix_to_pretty_latex(H)
    # Mapa zamiany nazw symboli na ładniejsze LaTeX-owe
    repl = Dict(
        :WntDiff => "W_D",
        :WntLoc  => "W_L",
        :DkkC    => "C",
        :DKKC    => "C",
        :β2      => "\\beta_{2}",
        :β3      => "\\beta_{3}",
        :β5      => "\\beta_{5}",
        :β6      => "\\beta_{6}",
    )

    # LaTeX samej macierzy w środowisku pmatrix
    pm = latexify(H;
        mathtype = :pmatrix,
        symbol_replace = repl
    )

    # Opakowanie w align*
    return "\\begin{align*}\n" * String(pm) * "\n\\end{align*}"
end

println(matrix_to_pretty_latex(H))