module Struktury

    export Parameters, Coefficients, Diffusions, VariablesVector

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

    export ApplyLaplacian, NonlinearityFunction

    # Different Variants of Nonlinearities

    #Variant 1
    function N1(Par::Parameters,Var::VariablesVector, t::Float64) 
        return VariablesVector(
                                Var.SD  .* Par.Coef.β6 ./ ((1 .+ Var.DkkA).*(1 .+ Var.DkkC).*(1 .+ Par.Coef.β3 .* Var.WntLoc)) .- Var.WntLoc,
                                Par.Coef.β1 ./ (1 .+ Par.Coef.β4 .* Var.WntLoc) .- Var.DkkA,
                                Par.Coef.β2 ./ 2.0 .* Var.WntLoc .+ Par.Coef.β2 ./ 2.0 .* Var.SD .- Var.WntDiff,
                                Var.WntDiff ./ (1 .+ Par.Coef.β5 .* Var.WntLoc) .- Var.DkkC,
                                Var.WntLoc .- Var.SD,
                              );
    end

    #Variant 2 with u^2 instead of exponents       
    function N2(Par::Parameters,Var::VariablesVector) 
        return nothing
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

    import Base: +


    +(x::T, y::T) where T = T([getfield(x,i) + getfield(y,i) for i in fieldnames(T)]...)



    ν = [0.0, 3.8154e-03, 0.4433, 6.0713e-03, 0.0004];
    β = [1.0629, 540.4003, 1.1596, 11.5964, 11.5964, 4.8254];

    U0 = [0.0026, 1.0317, 1.40926, 1.36789, 0.0026];

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
    Ini = IniCst   # Initial Conditions


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


