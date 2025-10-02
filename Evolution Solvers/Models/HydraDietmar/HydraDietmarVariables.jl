
using ..Struktury
using ..SimParam

dt = 0.01; # Time step
# dt = 1.0;
# dt = 0.1;

D = 0.001; # Diffusion Coefficient
κ = 10.0; # Nonlinearity Coefficient

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
Set = Parameters(Diffusions(D), Coefficients(κ)); # Parameters

# Set = Parameters(Diffusions(D), Coefficients(κ)); # Parameters
# VIni = rand(SimParam.N); # Initial Conditions
# VIni[091:100] = 10.0*ones(10)*1.0;
# VIni = VariablesVector(VIni); # Initial Conditions
