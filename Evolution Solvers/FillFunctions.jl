module FillMatrix

    using LinearAlgebra
    using ..SimParam

    export FiniteDiffercePeriodic, FiniteDifferceNeumann, FiniteDifferceDirichlet
    
    function FiniteDiffercePeriodic(Var::Vector{Float64})
        L = zeros(SimParam.N, SimParam.N)
        k =  floor(Int, length(Var)/2 - 0.5);
        for i = 1:SimParam.N
            for j = -k:1:k
                L[i,mod(i-j-1,SimParam.N)+1] = Var[k+1-j];
            end
        end

        return L
    end

    function FiniteDifferceNeumann(Var::Vector{Float64})
        L = zeros(SimParam.N, SimParam.N)
        k =  floor(Int, length(Var)/2 - 0.5);
        for i = 1:SimParam.N
            L[i,i] = Var[k+1];
            for j = 1:k
                L[i,i-j>0 ? i-j : 1 ] = Var[k+1-j];
                L[i,i+j <SimParam.N ? i+j : SimParam.N] = Var[k+1+j];
            end
        end
        for i=1:k
            L[i, 1] = sum(Var[1:k+2-i])
            L[end-i+1, end] = sum(Var[k+i:end])
        end
        return L
    end

    function FiniteDifferceDirichlet(Var::Vector{Float64})
        L = zeros(SimParam.N, SimParam.N)
        k =  floor(Int, length(Var)/2 - 0.5);
        for i = 1:SimParam.N
            for j = -min(k,i-1):1:0
                L[i+j,i] = Var[k+1+j];
            end
            for j = 0:1:min(k,SimParam.N-i)
                L[i+j,i] = Var[k+1+j];
            end
        end
        return L
    end 
    
    # Same functions but with N as argument

    function FiniteDiffercePeriodic(Var::Vector{Float64}, N::Int64)
        L = zeros(N, N)
        k =  floor(Int, length(Var)/2 - 0.5);
        for i = 1:N
            for j = -k:1:k
                L[i,mod(i-j-1,N)+1] = Var[k+1-j];
            end
        end

        return L
    end

    function FiniteDifferceNeumann(Var::Vector{Float64}, N::Int64)
        L = zeros(N, N)
        k =  floor(Int, length(Var)/2 - 0.5);
        for i = 1:N
            L[i,i] = Var[k+1];
            for j = 1:k
                L[i,i-j>0 ? i-j : 1 ] = Var[k+1-j];
                L[i,i+j <N ? i+j : N] = Var[k+1+j];
            end
        end
        for i=1:k
            L[i, 1] = sum(Var[1:k+2-i])
            L[end-i+1, end] = sum(Var[k+i:end])
        end
        return L
    end

    function FiniteDifferceDirichlet(Var::Vector{Float64}, N::Int64)
        L = zeros(N, N)
        k =  floor(Int, length(Var)/2 - 0.5);
        for i = 1:N
            for j = -min(k,i-1):1:0
                L[i+j,i] = Var[k+1+j];
            end
            for j = 0:1:min(k,N-i)
                L[i+j,i] = Var[k+1+j];
            end
        end
        return L
    end 

end