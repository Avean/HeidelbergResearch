
using ..Struktury
using ..SimParam

include("ReceptorBasedSuplementary.jl")

# dt = 0.000001; # Time step
dt = 0.5;
# dt = 1.0;

D0 = 0.0; # Diffusion Coefficient
D1 = 0.015; # Diffusion Coefficient
D2 = 0.1; # Diffusion Coefficient

u1 = 1.0; # Nonlinearity Coefficient
u2 = 1.0;
u3 = 0.6;

m1 = 2.5; # Nonlinearity Coefficient
m2 = 9.68;
m3 = 7.0;


NonlinearityFunction = "Nonlinearity 1"; # Classical model
# NonlinearityFunction = "Nonlinearity 2"; # Shadow limit
# NonlinearityFunction = "Nonlinearity 3"; # Fast reaction for V - one equation
# NonlinearityFunction = "Nonlinearity 4"; # Nonlinearity Variant with u^10
# NonlinearityFunction = "Nonlinearity 5"; # Nonlinearity Variant with u^2 / int u^3

# VIni = VariablesVector(u_p.*ones(SimParam.N) + 0.01.*rand(SimParam.N), 
#                        v_p.*ones(SimParam.N)+ 0.01.*rand(SimParam.N),
#                        w_p.*ones(SimParam.N)+ 0.01.*rand(SimParam.N)); # Initial Conditions


VIni = VariablesVector(0.0.*ones(SimParam.N) + 0.0.*rand(SimParam.N), 
                       0.0.*ones(SimParam.N)+ 0.0.*rand(SimParam.N),
                       0.0.*ones(SimParam.N)+ 0.0.*rand(SimParam.N)); # Initial Conditions


# Width = 60;

# VIni.u[1:Width] = u_p * ones(Width);
# VIni.v[1:Width] = v_p * ones(Width);
# VIni.w[1:Width] = w_p * ones(Width);

Location1 = 110;
Width1 = 640;

VIni.u[Location1:Location1 + Width1 - 1] = u_p * ones(Width1);
VIni.v[Location1:Location1 + Width1 - 1] = v_p * ones(Width1);
VIni.w[Location1:Location1 + Width1 - 1] = w_p * ones(Width1);


# Location1 = 200;
# Width1 = 70;

# VIni.u[Location1:Location1 + Width1 - 1] = u_p * ones(Width1);
# VIni.v[Location1:Location1 + Width1 - 1] = v_p * ones(Width1);
# VIni.w[Location1:Location1 + Width1 - 1] = w_p * ones(Width1);


# Location2 = 700;
# Width2 = 300;

# VIni.u[Location2:Location2 + Width2 - 1] = 0.0 * ones(Width2);
# VIni.v[Location2:Location2 + Width2 - 1] = 0.0 * ones(Width2);
# VIni.w[Location2:Location2 + Width2 - 1] = 0.0 * ones(Width2);

# Location3 = 390;
# Width3 = 30;

# VIni.u[Location3:Location3 + Width3 - 1] = 1.048 * ones(Width3);
# VIni.v[Location3:Location3 + Width3 - 1] = 0.68 * ones(Width3);
# VIni.w[Location3:Location3 + Width3 - 1] = 4.89 * ones(Width3);

# Location4 = 800;

# VIni.u[Location4:Location4 + Width - 1] = 1.048 * ones(Width);
# VIni.v[Location4:Location4 + Width - 1] = 0.68 * ones(Width);
# VIni.w[Location4:Location4 + Width - 1] = 4.89 * ones(Width);

using ..Struktury

Set = Parameters(Diffusions(D0, D1,D2), Coefficients(u1,u2,u3, m1, m2, m3)); # Parameters

