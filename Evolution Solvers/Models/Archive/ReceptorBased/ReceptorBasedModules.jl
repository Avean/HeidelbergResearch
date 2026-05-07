module Struktury

    export Parameters, Coefficients, Diffusions, VariablesVector

    struct Coefficients
        u1::Float64
        u2::Float64
        u3::Float64

        m1::Float64
        m2::Float64
        m3::Float64        
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
        # Add more variables as needed
    end


end

module SimParam
    N = 1000; # Number of discretization points
    dt = 0.5; # Time step
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
                                - Par.Coef.u1 .* Var.u .+ Par.Coef.m1 .* Var.u .* Var.v ./ (1 .+ Var.u.*Var.v), 
                                - Par.Coef.u2 .* Var.v .+ Par.Coef.m2 .* Var.u .* Var.v ./ (1 .+ Var.u.*Var.v) - Var.v .* Var.w,
                                - Par.Coef.u3 .* Var.w .+ Par.Coef.m3 .* Var.u .* Var.v ./ (1 .+ Var.u.*Var.v),
                              );
    end

    #Variant 2 with u^2 instead of exponents       
    function N2(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                                0.0,
                                0.0,
                                0.0
                              );
    end

    #Variant 3 with u^3 instead of exponents       
    function N3(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                                0.0,
                                0.0,
                                0.0
                                );
    end
    
    function N4(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                                0.0,
                                0.0,
                                0.0
                              );
    end
    
    function N5(Par::Parameters,Var::VariablesVector) 
        return VariablesVector(
                            0.0,
                            0.0,
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

    export  CreateTower
       
    SetL1 = Parameters(Diffusions(0.0,0.06,0.2), Coefficients(1.0, 1.0, 0.6, 2.5, 9.68, 7.0)); # Parameters
    SetL2 = Parameters(Diffusions(0.0,0.03,0.1), Coefficients(1.0, 1.0, 0.6, 2.5, 9.68, 7.0)); # Parameters
    SetL3 = Parameters(Diffusions(0.0,0.012,0.03), Coefficients(1.0, 1.0, 0.6, 2.5, 9.68, 7.0)); # Parameters
    SetL12 = Parameters(Diffusions(0.0,0.03,0.2), Coefficients(1.0, 1.0, 0.6, 2.5, 9.68, 7.0)); # Parameters
    SetL23 = Parameters(Diffusions(0.0,0.015,0.1), Coefficients(1.0, 1.0, 0.6, 2.5, 9.68, 7.0)); # Parameters
    SetL123 = Parameters(Diffusions(0.0,0.015,0.2), Coefficients(1.0, 1.0, 0.6, 2.5, 9.68, 7.0)); # Parameters                        
    
    include("ReceptorBasedSuplementary.jl")

    # Set2 = Parameters(Diffusions(0.2,0.2), Coefficients(10.0));
    # Set3 = Parameters(Diffusions(0.25,0.25), Coefficients(10.0));
    # Set4 = Parameters(Diffusions(0.3,0.3), Coefficients(10.0));
    # Set5 = Parameters(Diffusions(0.5,0.5), Coefficients(10.0));

    ### Initial contidion - perturben constant
    VIniCst = VariablesVector(u_p .*ones(SimParam.N) + 0.1 .*(rand(SimParam.N).-0.5), 
                              v_p .*ones(SimParam.N) + 0.1 .* (rand(SimParam.N).-0.5),
                              w_p .*ones(SimParam.N) + 0.1 .* (rand(SimParam.N).-0.5)
                              ); # Initial Conditions

    function CreateTower(Location::Vector{Vector{Float64}}, u::Vector{Float64}, u_p::Float64)
        for Loc in Location
            Loc = floor.(Int,Loc.*SimParam.N);
            u[Loc[1]:Loc[2]] = u_p .* ones(Loc[2]-Loc[1] + 1);
        end
        return u
    end

    ### Initial contidion - Single tower
    Location = [[0.2, 0.4]];

    u0 = zeros(SimParam.N);
    v0 = zeros(SimParam.N);
    w0 = zeros(SimParam.N);

    u0 = CreateTower(Location, u0, u_p);
    v0 = CreateTower(Location, v0, v_p);    
    w0 = CreateTower(Location, w0, w_p);
    
    VIniTower = VariablesVector(u0, 
                                v0,
                                w0
                                ); # Initial Conditions
    ### Initial contidion - Double tower
    Location = [
                [0.2, 0.4],
                [0.6, 0.8]
                ];

    u0 = zeros(SimParam.N);
    v0 = zeros(SimParam.N);
    w0 = zeros(SimParam.N);

    u0 = CreateTower(Location, u0, u_p);
    v0 = CreateTower(Location, v0, v_p);    
    w0 = CreateTower(Location, w0, w_p);
    
    VIniDoubleTower = VariablesVector(u0, 
                                v0,
                                w0
                                ); # Initial Conditions

    ### Initial contidion - Single tower perturbed
    Location = [[0.2, 0.4]];

    u0 = 0.01.*rand(SimParam.N);
    v0 = 0.01.*rand(SimParam.N);
    w0 = 0.01.*rand(SimParam.N);

    u0 = CreateTower(Location, u0, u_p);
    v0 = CreateTower(Location, v0, v_p);    
    w0 = CreateTower(Location, w0, w_p);
    
    VIniTowerPerturbed = VariablesVector(u0, 
                                v0,
                                w0
                                ); # Initial Conditions

    ### Initial contidion - Double tower perturbed                            
    Location = [
                [0.2, 0.4],
                [0.6, 0.8]
                ];

    u0 = 0.01.*rand(SimParam.N);
    v0 = 0.01.*rand(SimParam.N);
    w0 = 0.01.*rand(SimParam.N);

    u0 = CreateTower(Location, u0, u_p);
    v0 = CreateTower(Location, v0, v_p);    
    w0 = CreateTower(Location, w0, w_p);
    
    VIniDoubleTowerPerturbed = VariablesVector(u0, 
                                v0,
                                w0
                                ); # Initial Conditions
end


