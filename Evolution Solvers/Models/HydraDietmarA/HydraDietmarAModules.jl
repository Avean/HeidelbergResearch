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
        u::Vector{Float64}
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
                                - Var.u + (Par.Coef.κ  +  TimeSlope(Par,t)) * exp.(Var.u) ./ mean(exp.(Var.u)) 
                              );
    end


    #Variant 2 with u^2 instead of exponents       
    function N2(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                                - Var.u + Par.Coef.κ * (Var.u).^2 ./ mean((Var.u).^2)
                              );
    end

    #Variant 3 with u^3 instead of exponents       
    function N3(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                                - Var.u + Par.Coef.κ * (Var.u).^3 ./ mean((Var.u).^3)
                              );
    end
    
    function N4(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                                - Var.u + Par.Coef.κ * (Var.u).^10 ./ mean((Var.u).^10)
                              );
    end
    
    function N5(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                            - Var.u + Par.Coef.κ * (Var.u).^2 ./ mean((Var.u).^3)
                              );
    end


    ###### 
    ### Nonlinearities with kernels
    ######

    
    ### Simple Rectangle Kernel ####



    function N6(Par::Parameters,Var::VariablesVector) 
        KernelSize = 0.25;
        One = [ones(map(Int,SimParam.N*KernelSize)); 1; ones(map(Int,SimParam.N*KernelSize))]; 
        M = FiniteDiffercePeriodic(One)./sum(One);
        return VariablesVector(
                            - Var.u + Par.Coef.κ * exp.(Var.u) ./ (M*exp.(Var.u))
                              );
    end

    ### Cosine Kernel ####
    

    function N7(Par::Parameters, Var::VariablesVector)
        KernelSize = 1/5;
        C = range(0.0,pi/2,map(Int,(map(Int,SimParam.N*KernelSize))));
        CosKernel = [reverse(cos.(C)); 1; cos.(C)];
        CosKernel = CosKernel ./ sum(CosKernel);
        M_cos = FiniteDiffercePeriodic(CosKernel);
        return VariablesVector(
            - Var.u + Par.Coef.κ * exp.(Var.u) ./ (M_cos * exp.(Var.u))
            );
    end


    ### Gaussian Kernel ####
    
    function N8(Par::Parameters, Var::VariablesVector)
        KernelSize = 0.5;   
        C = range(0.0,4,map(Int,(map(Int,SimParam.N*KernelSize))));
        GaussKernel = [reverse(exp.(-C.^2)); 1; exp.(-C.^2)];
        GaussKernel = GaussKernel ./ sum(GaussKernel);
        M_gauss = FiniteDiffercePeriodic(GaussKernel);
        return VariablesVector(
                            - Var.u + Par.Coef.κ * exp.(Var.u) ./ (M_gauss * exp.(Var.u))
                              );
    end

    
    ### Introduce MORE kernels ####
    ### Kernel structure with matrix and parameter a ##

    ## Kernel only in the nominantor ## 

    function N9(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                            - Var.u + Par.Coef.κ * (Par.Ker.M*exp.(Var.u)) ./ mean(exp.(Var.u))
                              );
    end

    ## Kernel in both nominator and denominator with parameter a ##

    function N10(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                            - Var.u + Par.Coef.κ * ((Par.Ker.M*exp.(sign(Par.Ker.a).*Var.u)).^(Par.Ker.a)) ./ mean((Par.Ker.M*exp.(sign(Par.Ker.a).*Var.u)).^(Par.Ker.a))
                              );
    end

    


    NonlinearityFunction = Dict(
        "Nonlinearity 1" => N1,
        "Nonlinearity 2" => N2,
        "Nonlinearity 3" => N3,
        "Nonlinearity 4" => N4,
        "Nonlinearity 5" => N5,
        "Nonlinearity 6" => N6,
        "Nonlinearity 7" => N7,
        "Nonlinearity 8" => N8,
        "Nonlinearity 9" => N9,
        "Nonlinearity 10" => N10
        )
end


module Sets
    using ..Struktury
    using ..SimParam
    using ..Nonlinearity

    import Base: +



    +(x::T, y::T) where T = T([getfield(x,i) + getfield(y,i) for i in fieldnames(T)]...)

    

    D = [1e-4];
    κ = 0.0;

    
    ####### Find automatic Steady state ##########
    # U0, JacC,  = SetNonlinearity(Sets.NonlinearityType, Sets.β)
    JacC = [1.0];
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

    function VectorToIni(V::Float64)
        W = V .* ones(SimParam.N) 
        return VariablesVector(W)
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

    IniCst = VectorToIni(κ[1])   # Constant initial condition

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




