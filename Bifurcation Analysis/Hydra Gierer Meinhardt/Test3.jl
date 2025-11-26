using BenchmarkTools





X = rand(10000)

function  LapFor!(y,x)
    @inbounds @simd for k in 2:length(x)-1
        y[k] = x[k-1] + x[k+1] - 2*x[k]
    end
    return y
end


function LapVec(X)
    return  Y = X[1:end-2] + X[3:end] - 2*X[2:end-1] 
end

function LapTridiag(x,L)
    return L * x
end


Y = zeros(size(X));
    n = length(X)
    L = Tridiagonal(
        fill(1.0, n-1),   # poddiagonala
        fill(-2.0, n),    # główna
        fill(1.0, n-1)    # naddiagonala
    )
@btime LapVec(X);
@btime LapFor!(Y,X);
@btime LapTridiag(X,L);





using LinearAlgebra
using BenchmarkTools

 n = 10_000          # ważne: const, żeby kompilator dostał konkretny typ
 x = rand(n)
y = similar(x)

L = Tridiagonal(    # macierz budujemy raz, nie w każdym pomiarze
    ones(n-1),
    -2 .* ones(n),
    ones(n-1)
)

function lap_loop!(y, x)
    @assert length(y) == length(x)
    n = length(x)
    @inbounds @simd for i in 2:n-1
        y[i] = x[i-1] - 2x[i] + x[i+1]
    end
    # brzegi jak chcesz, np. Dirichlet 0:
    @inbounds begin
        y[1] = 0
        y[n] = 0
    end
    return y
end

@btime lap_loop!(Y, X);

function lap_vec!(y, x)
    @views y[2:end-1] .= x[1:end-2] .- 2 .* x[2:end-1] .+ x[3:end]
    @inbounds begin
        y[1] = 0
        y[end] = 0
    end
    return y
end

@btime lap_vec!($y, $x);

@btime $L * $x;





w = range(0,1,100)
x1 = sin.(w)
x2 = cos.(w)
x3 = exp.(w)
x4 = exp.(-w)
x5 = w.^2
t = 0;

x = [x1;x2;x3;x4;x5]

F = similar(x)

FGMeinhardt!(F,x,par_GM,t)

dx = 1
    for i in 2:n-1
        WntLOCb[i]  = ν1 * Lap2(WntLOC,i)/dx^2
        DkkAb[i]    = ν2 * Lap2(DkkA,i)/dx^2
        WntDiffb[i] = ν3 * Lap2(WntDiff,i)/dx^2;
        DkkCb[i]    = ν4 * Lap2(DkkC,i)/dx^2;
        SDb[i]      = ν5 * Lap2(SD,i)/dx^2;
    end
