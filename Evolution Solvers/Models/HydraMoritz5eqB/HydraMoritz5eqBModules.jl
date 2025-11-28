module Struktury

    export Parameters, Coefficients, Diffusions, VariablesVector, UnpackStruct

    mutable struct Coefficients
        β1::Float64
        β2::Float64
        β3::Float64
        β4::Float64
        β5::Float64
        β6::Float64
    end

    mutable struct Diffusions
        ν1::Float64
        ν2::Float64
        ν3::Float64
        ν4::Float64
        ν5::Float64
    end

    struct Parameters
        Diff::Diffusions
        Coef::Coefficients
    end

    struct VariablesVector
        WntLoc::Vector{Float64}
        DkkA::Vector{Float64}
        WntDiff::Vector{Float64}
        DkkC::Vector{Float64}
        SD::Vector{Float64}
    end

    function UnpackStruct(s)
        return [getfield(s, i) for i in fieldnames(typeof(s))]
    end
end

module SimParam
    N = 1000; # Number of discretization points
    dt = 0.1; # Time step
    L = 1; # Domain size
    dx = L/(N-1); # Spatial step
    x = range(0,L,N); # Discretisation Vector
    SicCosNodes = 20; # Sine Cosine Nodes
end




module  Nonlinearity
    using ..Struktury
    using ..SimParam

    using Statistics

    export ApplyLaplacian, NonlinearityFunction, SetNonlinearity, PrepareNonlinearity


    module SymbolicData

        using ..Struktury
        using ModelingToolkit
        using NLsolve


        @variables WntLoc DkkA WntDiff DkkC SD
        @parameters β1 β2 β3 β4 β5 β6

        Var = [WntLoc DkkA WntDiff DkkC SD];
        β = [β1 β2 β3 β4 β5 β6];
        Configuration = []
        X0 = ones(fieldcount(VariablesVector)) # Initial guess for constant solution


        function PreCalc(Fun::T) where T
            FunNonlinearity, = ModelingToolkit.build_function(Fun,Var, β, expression = Val(false))

            H =  ModelingToolkit.jacobian(Fun,Var);
            JacNonlinearity,  = ModelingToolkit.build_function(H, Var, β,expression = Val(false));
            
            return FunNonlinearity, JacNonlinearity
        end

        function EvaluateConstant(FunJac, β, X0 = SymbolicData.X0)
            display(X0)
            Fun, Jac = FunJac()
            Xs = nlsolve(x-> Fun(x, β), 
                        x-> Jac(x, β),
                        X0)
            return Xs.zero
        end

        function Init5eq() #Full Model
            
            ############# Nonlinearities ############
            F1 = β6 *SD / (1+DkkA) / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
            F2 = β1 / (1+ β4 * WntLoc) - DkkA
            F3 = β2 * WntLoc * SD - WntDiff
            F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
            F5 = WntLoc - SD            
            ########################################

            F = [F1 F2 F3 F4 F5]

            Configuration = [1, 1, 1, 1, 1]
            return PreCalc(F)
        end

        function Init5eqA() #5 equations, linear WntDiff
            
            ############# Nonlinearities ############
            F1 = β6 *SD / (1+DkkA) / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
            F2 = β1 / (1+ β4 * WntLoc) - DkkA
            F3 = β2 * WntLoc  - WntDiff
            F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
            F5 = WntLoc - SD            
            ########################################

            F = [F1 F2 F3 F4 F5]
            Configuration = [1, 1, 1, 1, 1]
            return PreCalc(F)
        end

        function Init4eqA() # Kick DkkA 
            
            ############# Nonlinearities ############
            F1 = β6 *SD  / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
            F2 = 0
            F3 = β2 * WntLoc * SD - WntDiff
            F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
            F5 = WntLoc - SD            
            ########################################

            F = [F1 F2 F3 F4 F5]
            Configuration = [1, 0, 1, 1, 1]
            return PreCalc(F)
        end

        function Init4eqB() # Kick DkkA and multiplication for WntDiff

            ############# Nonlinearities ############
            F1 = β6 *SD  / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
            F2 = 0
            F3 = β2 * WntLoc  - WntDiff
            F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
            F5 = WntLoc - SD            
            ########################################

            F = [F1 F2 F3 F4 F5]
            Configuration = [1, 0, 1, 1, 1]
            return PreCalc(F)
        end

        function Init4eqC() # Kick SD, quasi steady state

             ############# Nonlinearities ############
            F1 = β6 * WntLoc / (1+DkkA) / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
            F2 = β1 / (1+ β4 * WntLoc) - DkkA
            F3 = β2 * WntLoc * WntLoc - WntDiff
            F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
            F5 = 0           
            ########################################

            F = [F1 F2 F3 F4 F5]
            Configuration = [1, 1, 1, 1, 0]
            return PreCalc(F)
        end

        function Init4eqD() # Kick SD, WntDiff linear

             ############# Nonlinearities ############
            F1 = β6 * WntLoc / (1+DkkA) / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
            F2 = β1 / (1+ β4 * WntLoc) - DkkA
            F3 = β2 * WntLoc  - WntDiff
            F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
            F5 = 0           
            ########################################

            F = [F1 F2 F3 F4 F5]
            Configuration = [1, 1, 1, 1, 0]
            return PreCalc(F)
        end
        

        function Init4eqE() # Kick SD, WntDiff linear, substitue WntDiff in first equation

             ############# Nonlinearities ############
            F1 = β6 * WntDiff / (1+DkkA) / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
            F2 = β1 / (1+ β4 * WntLoc) - DkkA
            F3 = β2 * WntLoc  - WntDiff
            F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
            F5 = 0           
            ########################################

            F = [F1 F2 F3 F4 F5]
            Configuration = [1, 1, 1, 1, 0]
            return PreCalc(F)
        end

        function Init3eqA() # Kick SD, A and A linear WntDiff

            ############# Nonlinearities ############
            F1 = β6 * WntLoc  / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
            F2 = 0
            F3 = β2 * WntLoc  - WntDiff
            F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
            F5 = 0           
            ########################################

            F = [F1 F2 F3 F4 F5]
            Configuration = [1, 0, 1, 1, 0]
            return PreCalc(F)            
        end

        function Init3eqB() # Kick SD and A, linear WntDiff, substitute WntDiff for first equation

            ############# Nonlinearities ############
            F1 = β6 * WntDiff  / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
            F2 = 0
            F3 = β2 * WntLoc  - WntDiff
            F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
            F5 = 0           
            ########################################

            F = [F1 F2 F3 F4 F5]
            Configuration = [1, 0, 1, 1, 0]
            return PreCalc(F)
            
        end

        function Init3eqC() # Kick SD and, substitute WntDiff for first equation

            ############# Nonlinearities ############
            F1 = β6 * WntDiff / (1+DkkC) / (1+ β3 * WntLoc) - WntLoc
            F2 = 0
            F3 = β2 * WntLoc^2  - WntDiff
            F4 = WntDiff / (1 + β5 * WntLoc) - DkkC
            F5 = 0           
            ########################################

            F = [F1 F2 F3 F4 F5]
            Configuration = [1, 0, 1, 1, 0]
            # SymbolicData.X0 = [47.459090682404835, 0.0, 1217178.87756134, 2207.61535199080, 0.0]
            SymbolicData.X0 = [42.59036855933033, 0.0, 7981.33377368621, 16.2622428298969, 0.0]
            return PreCalc(F)
            
        end
    end

    ChoosenFunction =[] 

    function SetNonlinearity(Fun, β)
        U, J(U,β), Flag = PrepareNonlinearity(Fun, β)
        if Flag
            println("Positive constant solution found. Setting nonlinearity...")
            X, = Fun()
            Nonlinearity.ChoosenFunction = X
            return U, J(U,β)
        else 
            println("Can not find positive solution, returning zeros")
            return zeros(fieldcount(VariablesVector)), zeros(fieldcount(VariablesVector), fieldcount(VariablesVector))
        end   
    end

    function PrepareNonlinearity(Fun, β)
        F, J = Fun()
        Flag = false
        U = SymbolicData.EvaluateConstant(Fun,β)
        # for k in 1:1000
        #     X0 = 20.0.*rand(fieldcount(VariablesVector))
        #     # display(k)
        #     U = SymbolicData.EvaluateConstant(Fun,β, X0)
        #     if all(U[SymbolicData.Configuration .== 1] .> 0.0 )
        #         break
        #     end
        # end

        if all(U .> 0.0 )
            Flag = true
            return U, J(U,β), Flag
        else
            return [], [], false
        end
    end
    

    # Different Variants of Nonlinearities

    #Variant 1
    function N1(Par::Parameters,Var::VariablesVector, t::Float64) 
        return VariablesVector(
                                Var.SD  .* Par.Coef.β6 ./ ((1 .+ Var.DkkA).*(1 .+ Var.DkkC).*(1 .+ Par.Coef.β3 .* Var.WntLoc)) .- Var.WntLoc,
                                Par.Coef.β1 ./ (1 .+ Par.Coef.β4 .* Var.WntLoc) .- Var.DkkA,
                                Par.Coef.β2  .* Var.WntLoc .* Var.SD  .- Var.WntDiff,
                                Var.WntDiff ./ (1 .+ Par.Coef.β5 .* Var.WntLoc) .- Var.DkkC,
                                Var.WntLoc .- Var.SD,
                              );
    end

    #Variant 2 with u^2 instead of exponents       
    function N2(Par::Parameters,Var::VariablesVector, t::Float64)
        return VariablesVector(
                                Var.WntDiff  .* Par.Coef.β6 ./ ((1 .+ Var.DkkC).*(1 .+ Par.Coef.β3 .* Var.WntLoc)) .- Var.WntLoc,
                                .- Var.DkkA,
                                Par.Coef.β2  .* Var.WntLoc .* Var.WntLoc  .- Var.WntDiff,
                                Var.WntDiff ./ (1 .+ Par.Coef.β5 .* Var.WntLoc) .- Var.DkkC,
                                .- Var.SD,
                              );
    end

    #Variant 3 with u^3 instead of exponents       
    function N3(Par::Parameters,Var::VariablesVector) 
        return nothing
    end
    
    function N4(Par::Parameters,Var::VariablesVector) 
        return nothing
    end
    
    function N5(Par::Parameters,Var::VariablesVector) 
        return nothing
    end

    NonlinearityFunction = Dict(
        "Nonlinearity 1" => N1,
        "Nonlinearity 2" => N2,
        "Nonlinearity 3" => N3,
        "Nonlinearity 4" => N4,
        "Nonlinearity 5" => N5
    )
end


module Sets
    using ..Struktury
    using ..SimParam
    using ..Nonlinearity

    import Base: +



    +(x::T, y::T) where T = T([getfield(x,i) + getfield(y,i) for i in fieldnames(T)]...)

    ######### Choose nonlinearity here - last part #########
    NonlinearityType = Nonlinearity.SymbolicData.Init3eqC
    ########################################################

    ν = [0.0, 3.8154e-05, 0.4433, 6.0713e-08, 0.0004];
    # β = [1.0629, 540.4003, 1.1596, 11.5964, 11.5964, 4.8254];
    β =  [1.06, 4.4, 1.2, 11.5, 11.5, 4.8];
            
    U0 = []
    JacC = []
    X0 = ones(fieldcount(VariablesVector))

    Fun, Jac = NonlinearityType()

    ####### Find automatic Steady state ##########
    # U0, JacC,  = SetNonlinearity(Sets.NonlinearityType, Sets.β)
    ##############################################

    ####### Force Steady state ##########
    U0 = Nonlinearity.SymbolicData.X0
    JacC = Jac(U0, Sets.β)
    ##############################################


    XDDI = []

    function DisplayDDI()
        println("Unstable eigennodes of the linearized system")

        for i in eachindex(XDDI)
            println(XDDI[i]) 
        end
    end

    function VectorToIni(V::Vector{Float64})
        W = [v .* ones(SimParam.N) for v in V]
        return VariablesVector(W...)
    end

    function PerturbationRandom(l::Float64)
        W = [2.0 .* l.* (rand(SimParam.N).-0.5) for i in 1:fieldcount(VariablesVector)]
        return VariablesVector(W...)
    end 

    function PerturbationRandomPlus(l::Float64)
        W = [ l.* rand(SimParam.N) for i in 1:fieldcount(VariablesVector)]
        return VariablesVector(W...)
    end

    function SetInitialConditions(V::Vector)
        W = [V[i] for i in 1:fieldcount(VariablesVector)]
        return VariablesVector(W...)
    end

    IniCst = VectorToIni(U0)   # Constant initial condition

    IniCstPerturbed = [ # Initial conditions plus small pertubation around constant with amplitude (ℓ)
                        IniCst + PerturbationRandom(0.01),
                        IniCst + PerturbationRandom(0.03),
                        IniCst + PerturbationRandom(0.001)      
                      ]

    IniCstPerturbedPlus = [ # Initial conditions plus small positive perturbation with amplitude (ℓ) 
                        IniCst + PerturbationRandomPlus(0.1 ),
                        IniCst + PerturbationRandomPlus(0.3 ),
                        IniCst + PerturbationRandomPlus(0.7 ),
                        IniCst + PerturbationRandomPlus(1.0 ),
                        IniCst + PerturbationRandomPlus(3.0 ),
                        IniCst + PerturbationRandomPlus(10.0)
                      ]


    Par = Parameters(Diffusions(ν...), Coefficients(β...)); # Parameters

    Ini = IniCstPerturbed[3]; # Here Choose initial condition


    function ResetParameters()
        Sets.Par = Parameters(Diffusions(ν...), Coefficients(β...));
    end

    module BifurcationData

        # Diffusion D1
        Points1 = [-1.0, -2.0]
        Names1  = ["None", "None"]

        # Diffusion D2
        Points2 = [-1.0, -2.0]
        Names2  = ["None", "None"]

        # Diffusion D3
        Points3 =  [
                    2.24378,
                    1.06801
                ]


        Names3  =  [
                    "k1",
                    "k2"
                ]

        # Diffusion D4
        Points4 = [-1.0, -2.0]
        Names4  = ["None", "None"]

        # Diffusion D5
        Points5 = [-1.0, -2.0]
        Names5  = ["None", "None"]

        # Vectors
        PointsX = [Points1, Points2, Points3, Points4, Points5]
        Points0 = [zeros(length(x)) for x in PointsX]
        Names  = [Names1,  Names2,  Names3,  Names4,  Names5]

    end

    using .BifurcationData
end




