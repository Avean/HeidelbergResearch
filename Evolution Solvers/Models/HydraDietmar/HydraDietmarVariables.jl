
using ..Struktury
using ..SimParam

dt = 1.0; # Time step
# dt = 0.1;

D = 0.0001; # Diffusion Coefficient
κ = 15.0; # Nonlinearity Coefficient
# κ = 5.0;

# D = 1e-6; # Diffusion Coefficient
# κ = 0.7; # Nonlinear term

# κ = 2.0; # Nonlinear term
# D = 1/20/pi^2; # Diffusion Coefficient



# NonlinearityFunction = "Nonlinearity 1"; # Nonlinearity Variant with exp
# NonlinearityFunction = "Nonlinearity 2"; # Nonlinearity Variant with u^2
# NonlinearityFunction = "Nonlinearity 3"; # Nonlinearity Variant with u^3
# NonlinearityFunction = "Nonlinearity 4"; # Nonlinearity Variant with u^10
# NonlinearityFunction = "Nonlinearity 5"; # Nonlinearity Variant with u^2 / int u^3

# NonlinearityFunction = "Nonlinearity 6"; # Nonlinearity Variant with exp and kernel 1/2
# NonlinearityFunction = "Nonlinearity 7"; # Nonlinearity Variant with exp and kernel cos
# NonlinearityFunction = "Nonlinearity 8"; # Nonlinearity Variant with exp and kernel gauss
NonlinearityFunction = "Nonlinearity 9"; # MORE kernels


VIni = VariablesVector(20.0 .*rand(SimParam.N)+ ones(SimParam.N)*κ); # Initial Conditions

# Using nonlinearity N9:
# Exponent in E(u) = (f(u)*K)^a, where f(u) = exp(-sign(a)u):
a = 1;
# Kernel K in E(u) = (f(u)*K)^a:
include("../../FillFunctions.jl")
using .FillMatrix

# Simple rectangle
kernelsize = 0.35;
One = [ones(map(Int,SimParam.N*kernelsize)); 1; ones(map(Int,SimParam.N*kernelsize))]; 
M = FiniteDiffercePeriodic(One)./sum(One);
# Cos kernel
kernelsize_cos = 0.1;
C_cos = range(0.0,pi/2,map(Int,(map(Int,SimParam.N*kernelsize_cos))));
CosKernel = [reverse(cos.(C_cos)); 1; cos.(C_cos)];
CosKernel = CosKernel ./ sum(CosKernel);
M_cos = FiniteDiffercePeriodic(CosKernel);
# Gauß kernel
kernelsize_gauss = 0.05;
C_gauss = range(0.0,4,map(Int,(map(Int,SimParam.N*kernelsize_gauss))));
GaussKernel = [reverse(exp.(-C_gauss.^2)); 1; exp.(-C_gauss.^2)];
GaussKernel = GaussKernel ./ sum(GaussKernel);
M_gauss = FiniteDiffercePeriodic(GaussKernel);

Set = Parameters(Diffusions(D), Coefficients(κ), Exponent(a), Kernel(M_cos)); # Parameters

# Set = Parameters(Diffusions(D), Coefficients(κ)); # Parameters
# VIni = rand(SimParam.N); # Initial Conditions
# VIni[091:100] = 10.0*ones(10)*1.0;
# VIni = VariablesVector(VIni); # Initial Conditions
