module Struktury
    export Dyfuzje, ParametryB, ParametryC, ParametryK, Parametry, VariablesVector
    struct Dyfuzje
        a1::Float64  
        a2::Float64
        a3::Float64
        a4::Float64
        a5::Float64   
    end

    struct ParametryB
        b1::Float64  
        b2::Float64
        b3::Float64
        b4::Float64
        b5::Float64   
    end

    struct ParametryC
        c1::Float64  
        c2::Float64
        c3::Float64
        c4::Float64
        c5::Float64  
    end

    struct ParametryK
        k1::Float64  
        k2::Float64
        k3::Float64
        k4::Float64
        k5::Float64 
    end

    struct Parametry
        D::Dyfuzje
        B::ParametryB
        C::ParametryC
        K::ParametryK
    end

    struct VariablesVector
        WNTtoc::Vector{Float64}
        WNTpro::Vector{Float64}
        DKKa::Vector{Float64}
        DKKc::Vector{Float64}
        SD::Vector{Float64}
    end


end

module SimParam
    N = 1000; # Number of discretization points
    dt = 0.1; # Time step
    L = 1; # Domain size
    dx = L/N; # Spatial step
    x = range(0,L,N); # Discretisation Vector
    SicCosNodes = 3; # Sine Cosine Nodes
end


module Sets
    using ..Struktury
    using ..SimParam

    VIni1 = VariablesVector(ones(SimParam.N),ones(SimParam.N),ones(SimParam.N),ones(SimParam.N),ones(SimParam.N));

    D = Dyfuzje(1,1,1,1,1);
    PB = ParametryB(1,1,1,1,1);
    PC = ParametryC(1,1,1,1,1);
    PK = ParametryK(1,1,1,1,1);

    Set1 = Parametry(D,PB,PC,PK);

end


module  Nonlinearity
    using ..Struktury

    export NonlinearityFull,NonlinearityReduced1, ApplyLaplacian

    function NonlinearityFull(Par::Parametry,Var::VariablesVector) 
        a = Var.SD  .* Par.B.b1 ./ ((1 .+ Par.K.k1.*Var.DKKa).*(1 .+Var.DKKc.*Par.K.k2).*(1 .+Par.K.k3.*Var.WNTtoc)) .- Par.C.c1 .* Var.WNTtoc;
        b = Par.B.b2 ./ (1 .+Par.K.k4.*Var.WNTtoc) .- Par.C.c2.*Var.DKKa;
        c = Par.B.b3.*Var.SD .- Par.C.c3.*Var.WNTpro;
        d = Par.B.b4.*Var.WNTpro./(1 .+Par.K.k5.*Var.WNTtoc) .- Par.C.c4.*Var.DKKc;
        e = Par.B.b5.*Var.WNTtoc .- Par.C.c5.*Var.SD;
        return VariablesVector(a,b,c,d,e);
    end


    function NonlinearityReduced1(Par::Parametry,Var::VariablesVector) 
        a = Var.SD  * Par.B.b1 / ((1 + Par.K.k1*Var.DKKa)*(1+Var.DKKc*Par.K.k2)*(1+Par.K.k3*Var.WNTtoc)) - Par.C.c1 * Var.WNTtoc;
        b = Par.B.b2 / (1+Par.K.k4*Var.WNTtoc) - Par.C.c2*Var.DKKa;
        c = Par.B.b3*Var.SD - Par.C.c3*Var.WNTpro;
        d = Par.B.b4*Var.WNTpro/(1+Par.K.k5*Var.WNTtoc) - Par.C.c4*Var.DKKc;
        return VariablesVector(a,b,c,d,2);
    end

    function NonlinearityReduced2(a::Float64,b::Float64) 
        return a*b
    end

    function ApplyLaplacian(Var::VariablesVector, Lap)
        return VariablesVector(Lap(Var.WNTtoc), Lap(Var.WNTpro), Lap(Var.DKKa), Lap(Var.DKKa), Lap(Var.DKKc), Lap(Var.SD));
    end
end




