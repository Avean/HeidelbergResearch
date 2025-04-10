module Struktury

    export Parameters, Coefficients, Diffusions, VariablesVector

    struct Coefficients
        κ::Float64
    end

    struct Diffusions
        D1::Float64
        D2::Float64
    end

    struct Parameters
        Diff::Diffusions
        Coef::Coefficients
    end

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
    SicCosNodes = 200; # Sine Cosine Nodes
end

module  Nonlinearity
    using ..Struktury
    using  ..SimParam
    using Statistics

    export ApplyLaplacian, NonlinearityFunction

    # Different Variants of Nonlinearities

    #Variant 1
    function N1(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                                - Var.u + Par.Coef.κ * exp.(Var.u) ./ mean(exp.(Var.u)),
                                zeros(SimParam.N)
                              );
    end

    #Variant 2 with u^2 instead of exponents       
    function N2(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                                - Var.u + Par.Coef.κ * (Var.u).^2 ./ mean((Var.u).^2),
                                0.0
                              );
    end

    #Variant 3 with u^3 instead of exponents       
    function N3(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                                - Var.u + Par.Coef.κ * (Var.u).^3 ./ mean((Var.u).^3),
                                0.0
                              );
    end
    
    function N4(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                                - Var.u + Par.Coef.κ * (Var.u).^10 ./ mean((Var.u).^10),
                                0.0
                              );
    end
    
    function N5(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                            - Var.u + Par.Coef.κ * (Var.u).^2 ./ mean((Var.u).^3),
                            0.0
                              );
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

    VIni1 = VariablesVector(200 .*rand(SimParam.N),200 .*rand(SimParam.N));
    VIni2 = VariablesVector(1 .*rand(SimParam.N),1 .*rand(SimParam.N));
    VIni3 = VariablesVector(0.2 .*rand(SimParam.N).+10.0,0.000002 .*rand(SimParam.N).+10.0);

    Set1 = Parameters(Diffusions(0.1,0.1), Coefficients(10.0));
    Set2 = Parameters(Diffusions(0.2,0.2), Coefficients(10.0));
    Set3 = Parameters(Diffusions(0.25,0.25), Coefficients(10.0));
    Set4 = Parameters(Diffusions(0.3,0.3), Coefficients(10.0));
    Set5 = Parameters(Diffusions(0.5,0.5), Coefficients(10.0));
end