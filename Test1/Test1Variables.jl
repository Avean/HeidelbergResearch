
using ..Struktury
using ..SimParam

# dt = 0.0001; # Time step
dt = 1.0;
# dt = 0.01;

D = 0.01; # Diffusion Coefficient
κ = 10.0; # Nonlinearity Coefficient

# D = 1e-4; # Diffusion Coefficient
# κ = 0.5; # Nonlinear term

NonlinearityFunction = "Nonlinearity 1"; # Nonlinearity Variant with exp
# NonlinearityFunction = "Nonlinearity 2"; # Nonlinearity Variant with u^2
# NonlinearityFunction = "Nonlinearity 3"; # Nonlinearity Variant with u^3
# NonlinearityFunction = "Nonlinearity 4"; # Nonlinearity Variant with u^10
# NonlinearityFunction = "Nonlinearity 5"; # Nonlinearity Variant with u^2 / int u^3

VIni = VariablesVector(
                        0.002 .*rand(SimParam.N)+ ones(SimParam.N)*κ,
                        0.002 .*rand(SimParam.N)+ ones(SimParam.N)*κ
                        ); # Initial Conditions

Set = Parameters(Diffusions(D,D), Coefficients(κ)); # Parameters

