
using ..Struktury
using ..SimParam

# dt = 0.000001; # Time step
# dt = 1.0;
dt = 0.001;

D1 = 0.0001; # Diffusion Coefficient
D2 = 0.28; # Diffusion Coefficient

a = 2.0; # Nonlinearity Coefficient
b = 0.1;
c = 0.2; 

# NonlinearityFunction = "Nonlinearity 1"; # Classical model
# NonlinearityFunction = "Nonlinearity 2"; # Shadow limit
# NonlinearityFunction = "Nonlinearity 3"; # Fast reaction for V - one equation
# NonlinearityFunction = "Nonlinearity 4"; # Nonlinearity Variant with u^10
NonlinearityFunction = "Nonlinearity 5"; # Nonlinearity Variant with u^2 / int u^3

# VIni = VariablesVector(1.77 .*ones(SimParam.N) + 0.1 .*(rand(SimParam.N).-0.5), 
#                        15.77 .*ones(SimParam.N) + 0.1 .* (rand(SimParam.N).-0.5)
#                        ); # Initial Conditions

# Set = Parameters(Diffusions(D1,D2), Coefficients(a,b,c)); # Parameters

