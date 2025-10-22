includet("FillFunctions.jl")
using .FillMatrix
using ..Struktury
using ..SimParam

function KernelRectangle(KernelSize::Float64, a::Float64)
        One = [ones(map(Int,SimParam.N*KernelSize)); 1; ones(map(Int,SimParam.N*KernelSize))]; 
        M = FiniteDiffercePeriodic(One)./sum(One);
        return Kernels(M, a)
end


function KernelGaussian(Size::Float64, a::Float64)   
        C = range(0.0,1.0./Size*5,map(Int,(map(Int,SimParam.N*0.5))));
        GaussKernel = [reverse(exp.(-C.^2)); 1; exp.(-C.^2)];
        GaussKernel = GaussKernel ./ sum(GaussKernel);
        M = FiniteDiffercePeriodic(GaussKernel);
        return Kernels(M, a)
end


