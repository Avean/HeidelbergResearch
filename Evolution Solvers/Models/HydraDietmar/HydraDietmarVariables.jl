includet("../../Kernels.jl")
using ..Struktury
using ..SimParam

dt = 0.01; # Time step
# dt = 1.0;
# dt = 0.1;

# D = 0.01; # Diffusion Coefficient
# κ = 1.0; # Nonlinearity Coefficient

D = 1e-4; # Diffusion Coefficient
κ = 0.0; # Nonlinear term
Slope = 0.1; #Slope value
LCritical = 2.05;

# κ = 2.0; # Nonlinear term
# D = 1/20/pi^2; # Diffusion Coefficient



NonlinearityFunction = "Nonlinearity 1"; # Nonlinearity Variant with exp
# NonlinearityFunction = "Nonlinearity 2"; # Nonlinearity Variant with u^2
# NonlinearityFunction = "Nonlinearity 3"; # Nonlinearity Variant with u^3
# NonlinearityFunction = "Nonlinearity 4"; # Nonlinearity Variant with u^10
# NonlinearityFunction = "Nonlinearity 5"; # Nonlinearity Variant with u^2 / int u^3

# NonlinearityFunction = "Nonlinearity 6"; # Nonlinearity Variant with exp and kernel 1/2
# NonlinearityFunction = "Nonlinearity 7"; # Nonlinearity Variant with exp and kernel cos
# NonlinearityFunction = "Nonlinearity 8"; # Nonlinearity Variant with exp and kernel gauss
# NonlinearityFunction = "Nonlinearity 9"; # Kernel i the nominantor only
# NonlinearityFunction = "Nonlinearity 10"; # Kernel in both nominantor and denominator with parameter a

# Kernel = KernelRectangle(0.05, 1.0); #Rectangle Kernel (Size, a)
Kernel = KernelGaussian(0.5 , -5.0); #Gaussian Kernel (Size, a)
# plot(Kernel.M[500,:])

VIni = VariablesVector(2.0e-2 .*rand(SimParam.N)+ ones(SimParam.N)*κ, zeros(SimParam.N)); # Initial Conditions
Set = Parameters(Diffusions(D,0.0), Coefficients(κ,Slope,LCritical),Kernel); # Parameters


# Set = Parameters(Diffusions(D), Coefficients(κ)); # Parameters
# VIni = rand(SimParam.N); # Initial Conditions
# VIni[091:100] = 10.0*ones(10)*1.0;
# VIni = VariablesVector(VIni); # Initial Conditions

