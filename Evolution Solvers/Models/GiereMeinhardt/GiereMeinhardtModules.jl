module Struktury

    export Parameters, Coefficients, Diffusions, VariablesVector

    struct Coefficients
        a::Float64
        b::Float64
        c::Float64
    end

    struct Diffusions
        D1::Float64
        D2::Float64
        D3::Float64
    end

    struct Parameters
        Diff::Diffusions
        Coef::Coefficients
    end

    struct VariablesVector
        u::Vector{Float64}
        v::Vector{Float64}
        w::Vector{Float64}
    end


end

module SimParam
    N = 100; # Number of discretization points
    dt = 0.1; # Time step
    L = 1; # Domain size
    dx = L/N; # Spatial step
    x = range(0,L,N); # Discretisation Vector
    SicCosNodes = 20; # Sine Cosine Nodes
end




module  Nonlinearity
    using ..Struktury
    using  ..SimParam
    using Statistics

    export ApplyLaplacian, NonlinearityFunction

    # Different Variants of Nonlinearities

    # #Variant 1
    # function N1(Par::Parameters,Var::VariablesVector) 
    #     return VariablesVector(
    #                             Var.u.^2 ./ (Var.v .+ Par.Coef.a) - Par.Coef.b .* Var.u, 
    #                             Var.u.^2 - Par.Coef.c * Var.v
    #                           );
    # end

    # #Variant 2 with u^2 instead of exponents       
    # function N2(Par::Parameters,Var::VariablesVector) 
    #     return VariablesVector(
    #                             Var.u.^2 ./ (Var.v .+ Par.Coef.a) - Par.Coef.b .* Var.u, 
    #                             mean(Var.u.^2) .- Par.Coef.c * Var.v
    #                           );
    # end

    # #Variant 3 with u^3 instead of exponents       
    # function N3(Par::Parameters,Var::VariablesVector) 
    #     return VariablesVector(
    #                             Var.u.^2 ./ (mean(Var.u.^2) ./ Par.Coef.c .+ Par.Coef.a) - Par.Coef.b .* Var.u, 
    #                             zeros(SimParam.N)
    #                             );
    # end
    
    # function N4(Par::Parameters,Var::VariablesVector) 
    #     return VariablesVector(
    #                             0.0,
    #                             0.0
    #                           );
    # end
    
    function N5(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                            Var.v - Var.u,
                            Par.Coef.a .- (Par.Coef.b+2)*Var.v + Var.v.^2 .* Var.w + Var.u,
                            Par.Coef.b*Var.v - Var.v.^2 .* Var.w
                              );
    end

    NonlinearityFunction = Dict(
        # "Nonlinearity 1" => N1,
        # "Nonlinearity 2" => N2,
        # "Nonlinearity 3" => N3,
        # "Nonlinearity 4" => N4,
        "Nonlinearity 5" => N5
    )
end


module Sets
    using ..Struktury
    using ..SimParam

    # VIni1 = VariablesVector(1.77 .*ones(SimParam.N) + 0.1 .*(rand(SimParam.N).-0.5), 
    #                         15.77 .*ones(SimParam.N) + 0.1 .* (rand(SimParam.N).-0.5)
    #                         ); # Initial Conditions

    # VIni2 = VariablesVector(1.77 .*ones(SimParam.N) + 0.1 .*(rand(SimParam.N).-0.5), 
    #                         15.77 .*ones(SimParam.N)
    #                         ); # Initial Conditions


    # VIni3 = VariablesVector(1.0 .*rand(SimParam.N),
    #                         1.0 .*rand(SimParam.N));

    # VIni4 = VariablesVector(1.77.*ones(SimParam.N),
    #                         15.7.*ones(SimParam.N));

    p = 1.5;
    q = 3.5;

    VIni5 = VariablesVector(p.*ones(SimParam.N),
                            p.*ones(SimParam.N),
                            q/p.*ones(SimParam.N));

    Set1 = Parameters(Diffusions(0.0,0.3,1.0), Coefficients(p,q,1.0)); # Parameters

    # Set2 = Parameters(Diffusions(0.0001,0.002), Coefficients(2.0,0.1,0.2)); # Parameters


    # Set2 = Parameters(Diffusions(0.2,0.2), Coefficients(10.0));
    # Set3 = Parameters(Diffusions(0.25,0.25), Coefficients(10.0));
    # Set4 = Parameters(Diffusions(0.3,0.3), Coefficients(10.0));
    # Set5 = Parameters(Diffusions(0.5,0.5), Coefficients(10.0));
end


