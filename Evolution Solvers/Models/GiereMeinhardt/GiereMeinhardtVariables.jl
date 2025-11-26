using ..Struktury
using ..SimParam


######
###### Choose one of the following schemes
######

# Scheme = "ExpliciteEuler"; 
Scheme = "IMEX"; 
# Scheme = "SpectralSinCos"; 

######
###### Choose one of the following boundary conditions
######

# BoundaryConditions = "Dirichlet";
BoundaryConditions = "Neumann";
# BoundaryConditions = "Periodic";

######
###### Choose one of the following levels of Laplacian discretization
######

# Order = "2";
# Order = "4";
# Order = "6";
Order = "8";

#######
####### Setting Parameters
#######

# dt = 0.000001; # Time step
# dt = 1.0;
dt = 0.001;

D1 = 0.0001; # Diffusion Coefficient
D2 = 0.28; # Diffusion Coefficient

a = 2.0; # Nonlinearity Coefficient
b = 0.1;
c = 0.2; 

######
###### Choosing Nonlinearity
######

NonlinearityFunction = "Nonlinearity 1"; # Classical model

######
###### Initial conditions and parameters
######


# VIni = VariablesVector(1.77 .*ones(SimParam.N) + 0.1 .*(rand(SimParam.N).-0.5), 
#                        15.77 .*ones(SimParam.N) + 0.1 .* (rand(SimParam.N).-0.5)
#                        ); # Initial Conditions

Set = Parameters(Diffusions(D1,D2), Coefficients(a,b,c)); # Parameters




InitialConditions = VIni;
# InitialConditions = Sets.VIniTower;
# InitialConditions = Sets.CstStableTower(5.1, [0.5, 0.51]);
# InitialConditions = Sets.CstStableTowerRandom(5.0, [0.0, 1.0]);
# InitialConditions = Sets.CstStableMediumCstPerturb


ParameterSet = Set;
# ParameterSet = Sets.SetL2;
# ParameterSet = Sets.CstStable;

