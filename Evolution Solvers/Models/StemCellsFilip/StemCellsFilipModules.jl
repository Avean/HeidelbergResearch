module Struktury

    export Parameters, Coefficients, Diffusions, VariablesVector, Kernels, UnpackStruct

    Base.@kwdef mutable struct Coefficients
        r0::Float64
        b0::Float64
        p0::Float64
        β::Float64
        H::Float64
        K::Float64
    end

    mutable struct Diffusions
        DA::Float64
        DQ::Float64
    end

    struct Kernels
        M::Matrix{Float64}
        a::Float64
    end

    struct Parameters
        Diff::Diffusions
        Coef::Coefficients
    end

    struct VariablesVector
        Q::Vector{Float64}
        A::Vector{Float64}
    end


    function UnpackStruct(s)
        return [getfield(s, i) for i in fieldnames(typeof(s))]
    end
end

module SimParam
    N = 100; # Number of discretization points
    dt = 1e-1; # Time step
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

    function TimeSlope(Par::Parameters,t)
        return (Par.Coef.Slope * (t - Par.Coef.lbreak/Par.Coef.Slope * floor.(t*Par.Coef.Slope/Par.Coef.lbreak))).^1;
    end

    function r(Par::Parameters, Var::VariablesVector)
        return Par.Coef.r0 .* (Var.Q + Var.A) ./ (Par.Coef.K .+ Var.Q + Var.A)
    end

    function b(Par::Parameters, Var::VariablesVector)
        return Par.Coef.b0 ./ (1.0 .+ Par.Coef.β .* Var.Q)
    end

    function p(Par::Parameters, Var::VariablesVector)
        return Par.Coef.p0 .* (Var.A .+ Var.Q) ./ (Par.Coef.H .+ Var.Q + Var.A)
    end

    # Different Variants of Nonlinearities

    #Variant 1
    function N1(Par::Parameters,Var::VariablesVector,t::Float64) 
        return VariablesVector(
                                - r(Par,Var) .* Var.Q + 2.0 .* b(Par,Var) .* p(Par,Var) .* Var.A,
                                r(Par,Var) .* Var.Q - p(Par,Var) .* Var.A
                              );
    end

    #Variant 1
    function N2(Par::Parameters,Var::VariablesVector,t::Float64) 
        return VariablesVector(
                                - r(Par,Var) .* Var.Q + 2.0 .* b(Par,Var) .* p(Par,Var) .* Var.A,
                                r(Par,Var) .* Var.Q - p(Par,Var) .* Var.A
                              );
    end
    

    function N3(Par::Parameters,Var::VariablesVector,t::Float64) 
        return VariablesVector(
                                - r(Par,Var) .* Var.Q + 2.0 .* b(Par,Var) .* p(Par,Var) .* Var.A,
                                r(Par,Var) .* Var.Q - p(Par,Var) .* Var.A
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
    
    # r0 = 0.5; # Proliferation rate
    # b0 = 0.5; # Self-renewal probability
    # p0 = 0.5; # Differentiation rate
    # β = 0.5; # Feedback strength
    # H = 0.5; # Half-saturation constant
    # K = 0.5; # Carrying capacity

    r0 = 1.1462; # Proliferation rate
    b0 = 0.52446; # Self-renewal probability
    p0 = 1.99967; # Differentiation rate
    β = 0.00052702; # Feedback strength
    H = 246.39; # Half-saturation constant
    K = 66.808; # Carrying capacity

    Qcst = 100.0;
    Acst = 60.0;


    
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

    IniCst = VectorToIni([Qcst, Acst])   # Constant initial condition

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


    Par = Parameters(Diffusions(D...), Coefficients(r0, b0, p0, β, H, K)); # Parameters

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
        Sets.Par = Parameters(Diffusions(D...), Coefficients(r0, b0, p0, β, H, K));
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




