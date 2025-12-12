module Struktury

    export Parameters, Coefficients, Diffusions, VariablesVector, Kernels, UnpackStruct

    Base.@kwdef mutable struct Coefficients
        κ::Float64
        Slope::Float64
        lbreak::Float64
    end

    Coefficients(κ::Float64) = Coefficients(κ, 0.1, 2.1)

    mutable struct Diffusions
        D1::Float64
        D2::Float64
    end

    struct Kernels
        M::Matrix{Float64}
        a::Float64
    end

    struct Parameters
        Diff::Diffusions
        Coef::Coefficients
        Ker::Union{Kernels,Missing} 
    end

    struct VariablesVector
        WD::Vector{Float64}
        BK::Vector{Float64}
    end

    Parameters(D::Diffusions, C::Coefficients) = Parameters(D, C, missing)

    function UnpackStruct(s)
        return [getfield(s, i) for i in fieldnames(typeof(s))]
    end
end

module SimParam
    N = 1000; # Number of discretization points
    dt = 1e-2; # Time step
    L = 1; # Domain size
    dx = L/(N-1); # Spatial step
    x = range(0,L,N); # Discretisation Vector
    SicCosNodes = 100; # Sine Cosine Nodes
end




module  Nonlinearity
    using ..Struktury
    using ..SimParam

    using Statistics

    export ApplyLaplacian, NonlinearityFunction, SetNonlinearity, PrepareNonlinearity

function TimeSlope(Par::Parameters,t)
        return (Par.Coef.Slope * (t - Par.Coef.lbreak/Par.Coef.Slope * floor.(t*Par.Coef.Slope/Par.Coef.lbreak))).^1;
    end

    # Different Variants of Nonlinearities

    #Variant 1
    function N1(Par::Parameters,Var::VariablesVector,t::Float64) 
        return VariablesVector(
                                - Var.WD + Par.Coef.κ  * exp.(Var.BK) ./ mean(exp.(Var.BK)),
                                - Var.BK + Var.WD
                              );
    end

    #Variant 1
    function N2(Par::Parameters,Var::VariablesVector,t::Float64) 
        return VariablesVector(
                                - Var.WD + Par.Coef.κ  * exp.(Var.BK.^2) ./ mean(exp.(Var.BK.^2)),
                                - Var.BK + Var.WD
                              );
    end
    

    function N3(Par::Parameters,Var::VariablesVector,t::Float64) 
        return VariablesVector(
                                - Var.WD + 5.0 .* Var.BK .^2 - Var.WD, 
                                - Var.BK .+ 1.0 .* exp.(Var.BK) ./ mean(exp.(Var.WD)) ./ (1.0 .+ Var.BK)
                              );
    end

        function N4(Par::Parameters,Var::VariablesVector,t::Float64) 
        return VariablesVector(
                                - Var.WD + (Par.Coef.κ .+ Var.BK) .*  exp.(Var.WD) ./ mean(exp.(Var.WD)), 
                                - Var.BK .* 100.0 + 1.0 .* Var.WD
                              );
    end


    NonlinearityFunction = Dict(
        "Nonlinearity 1" => N1,
        "Nonlinearity 2" => N2,
        "Nonlinearity 3" => N3,
        )
end


module Sets
    using ..Struktury
    using ..SimParam
    using ..Nonlinearity

    import Base: +

    export SetTower


    +(x::T, y::T) where T = T([getfield(x,i) + getfield(y,i) for i in fieldnames(T)]...)

    

    D = [1e-4, 0.0]; # Diffusion Coefficients
    κ = 2.0;

    
    ####### Find automatic Steady state ##########
    # U0, JacC,  = SetNonlinearity(Sets.NonlinearityType, Sets.β)
    JacC = ones(fieldcount(VariablesVector), fieldcount(VariablesVector)) # Placeholder
    ##############################################


    XDDI = []
    Configuration = [1, 2];

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

    IniCst = VectorToIni([κ[1]/2, κ[1]/2])   # Constant initial condition

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


    Par = Parameters(Diffusions(D...), Coefficients(κ, 0.1, 2.2)); # Parameters

    Ini = IniCstPerturbed[3]; # Here Choose initial condition

    Perturbation = PerturbationRandom(0.0);

    P1 = PerturbationRandom(0.05);
    P2 = PerturbationRandom(0.1);
    P3 = PerturbationRandom(0.2);
    P4 = PerturbationRandom(0.5);
    P5 = PerturbationRandom(1.0);
    P6 = PerturbationRandom(2.0);
    P7 = PerturbationRandom(5.0);
    P8 = PerturbationRandom(10.0);
    
    function SetTower(x1, x2, amp)
        W = [zeros(SimParam.N) for i in 1:fieldcount(VariablesVector)]
        x1 = round(Int, x1 / SimParam.dx)
        x2 = round(Int, x2 / SimParam.dx)
        for i in 1:fieldcount(VariablesVector)
            W[i][x1:x2] .= amp
        end
        Sets.Perturbation =  VariablesVector(W...)
    end




    function PerturbState(U::VariablesVector)
        return U += Perturbation
    end


    function ResetParameters()
        Sets.Par = Parameters(Diffusions(D...), Coefficients(κ, 0.1, 2.2));
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




