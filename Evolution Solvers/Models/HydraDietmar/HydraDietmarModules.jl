module Struktury

    export Parameters, Coefficients, Diffusions, VariablesVector, Kernels

    struct Coefficients
        κ::Float64
    end

    struct Diffusions
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

    Parameters(D::Diffusions, C::Coefficients) = Parameters(D, C, missing)

    struct VariablesVector
        u::Vector{Float64}
        v::Vector{Float64}
        
    end
end

module SimParam
    N = 1000; # Number of discretization points
    L = 1; # Domain size
    dx = L/N; # Spatial step
    x = range(0,L,N); # Discretisation Vector
    SicCosNodes = 4; # Sine Cosine Nodes
end

module  Nonlinearity

    using  ..SimParam
    using ..Struktury
    using Statistics

    include("../../FillFunctions.jl")
    using .FillMatrix
    



    export ApplyLaplacian, NonlinearityFunction


    # Different Variants of Nonlinearities

    #Variant 1
    function N1(Par::Parameters,Var::VariablesVector,t::Float64) 
        return VariablesVector(
                                - Var.u + (Par.Coef.κ  +  5.05*(t/10 - floor.(t/10))) * exp.(Var.u) ./ mean(exp.(Var.u)),
                                   zeros(size(Var.u)) 
                                # - Var.u + (Par.Coef.κ  +  mean(Var.v).*1.0*(t/10 - floor.(t/10))) * exp.(Var.u) ./ mean(exp.(Var.u)),
                                # 0.2.*ones(SimParam.N).*(t/10 - floor.(t/10)) .- 0.03.*Var.v
                              );
    end


    #Variant 2 with u^2 instead of exponents       
    function N2(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                                - Var.u + Par.Coef.κ * (Var.u).^2 ./ mean((Var.u).^2),
                                zeros(size(Var.u))
                              );
    end

    #Variant 3 with u^3 instead of exponents       
    function N3(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                                - Var.u + Par.Coef.κ * (Var.u).^3 ./ mean((Var.u).^3),
                                zeros(size(Var.u))
                              );
    end
    
    function N4(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                                - Var.u + Par.Coef.κ * (Var.u).^10 ./ mean((Var.u).^10),
                                zeros(size(Var.u))
                              );
    end
    
    function N5(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                            - Var.u + Par.Coef.κ * (Var.u).^2 ./ mean((Var.u).^3),
                            zeros(size(Var.u))
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
                            - Var.u + Par.Coef.κ * exp.(Var.u) ./ (M*exp.(Var.u)),
                            zeros(size(Var.u))
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
            - Var.u + Par.Coef.κ * exp.(Var.u) ./ (M_cos * exp.(Var.u)),
            zeros(size(Var.u))
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
                            - Var.u + Par.Coef.κ * exp.(Var.u) ./ (M_gauss * exp.(Var.u)),
                            zeros(size(Var.u))
                              );
    end

    
    ### Introduce MORE kernels ####
    ### Kernel structure with matrix and parameter a ##

    ## Kernel only in the nominantor ## 

    function N9(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                            - Var.u + Par.Coef.κ * (Par.Ker.M*exp.(Var.u)) ./ mean(exp.(Var.u)),
                            zeros(size(Var.u))
                              );
    end

    ## Kernel in both nominator and denominator with parameter a ##

    function N10(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                            - Var.u + Par.Coef.κ * ((Par.Ker.M*exp.(sign(Par.Ker.a).*Var.u)).^(Par.Ker.a)) ./ mean((Par.Ker.M*exp.(sign(Par.Ker.a).*Var.u)).^(Par.Ker.a)),
                            zeros(size(Var.u))
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

    VIni1 = VariablesVector(200 .*rand(SimParam.N), zeros(SimParam.N));
    VIni2 = VariablesVector(1 .*rand(SimParam.N), zeros(SimParam.N));
    VIni3 = VariablesVector(0.000002 .*rand(SimParam.N).+10.0, zeros(SimParam.N));

    Set1 = Parameters(Diffusions(0.1,0.0), Coefficients(10.0));
    Set2 = Parameters(Diffusions(0.2,0.0), Coefficients(10.0));
    Set3 = Parameters(Diffusions(0.25,0.0), Coefficients(10.0));
    Set4 = Parameters(Diffusions(0.3,0.0), Coefficients(10.0));
    Set5 = Parameters(Diffusions(0.5,0.0), Coefficients(10.0));

    CstUnstable = Parameters(Diffusions(0.01,0.0), Coefficients(10.0));
    CstStable = Parameters(Diffusions(1e-6,0.0), Coefficients(0.7));

    CstStableSmallCstPerturb = VariablesVector(0.2 .*rand(SimParam.N) .+ CstStable.Coef.κ,zeros(SimParam.N));
    CstUnstableSmallCstPerturb = VariablesVector(0.2 .*rand(SimParam.N) .+ CstUnstable.Coef.κ,zeros(SimParam.N));

    CstStableMediumCstPerturb = VariablesVector(1.0 .*rand(SimParam.N) .+ CstStable.Coef.κ,zeros(SimParam.N));
    CstUnstableMediumCstPerturb = VariablesVector(1.0 .*rand(SimParam.N) .+ CstUnstable.Coef.κ,zeros(SimParam.N));   

    CstStableLargeCstPerturb = VariablesVector(20.0 .*rand(SimParam.N) .+ CstStable.Coef.κ,zeros(SimParam.N));
    CstUnstableLargeCstPerturb = VariablesVector(20.0 .*rand(SimParam.N) .+ CstUnstable.Coef.κ, zeros(SimParam.N));

    function CstStableTower(Height::Float64,Location::Vector{Float64})
        W = CstStable.Coef.κ .* ones(SimParam.N);
        Loc = floor.(Int,Location.*SimParam.N);
        Loc[1] = max(1,Loc[1]);
        Loc[2] = min(SimParam.N,Loc[2]);
        W[Loc[1]:Loc[2]] += Height .* ones(Loc[2]-Loc[1]+1);
        return VariablesVector(W,zeros(SimParam.N));
    end
    
    function CstUnstableTower(Height::Float64,Location::Vector{Float64})
        W = CstUnstable.Coef.κ .* ones(SimParam.N);
        Loc = floor.(Int,Location.*SimParam.N);
        Loc[1] = max(1,Loc[1]);
        Loc[2] = min(SimParam.N,Loc[2]);
        W[Loc[1]:Loc[2]] += Height .* ones(Loc[2]-Loc[1]+1);
        return VariablesVector(W,zeros(SimParam.N));
    end

    function CstStableTowerRandom(Height::Float64,Location::Vector{Float64})
        W = CstStable.Coef.κ .* ones(SimParam.N);
        Loc = floor.(Int,Location.*SimParam.N);
        Loc[1] = max(1,Loc[1]);
        Loc[2] = min(SimParam.N,Loc[2]);
        W[Loc[1]:Loc[2]] += Height .* rand(Loc[2]-Loc[1]+1);
        return VariablesVector(W,zeros(SimParam.N));
    end
    
    function CstUnstableTowerRandom(Height::Float64,Location::Vector{Float64})
        W = CstUnstable.Coef.κ .* ones(SimParam.N);
        Loc = floor.(Int,Location.*SimParam.N);
        Loc[1] = max(1,Loc[1]);
        Loc[2] = min(SimParam.N,Loc[2]);
        W[Loc[1]:Loc[2]] += Height .* rand(Loc[2]-Loc[1]+1);
        return VariablesVector(W,zeros(SimParam.N));
    end
   
end