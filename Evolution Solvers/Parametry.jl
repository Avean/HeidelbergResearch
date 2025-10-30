
include("Models/"*ModelName*"/"*ModelName*"Variables.jl")

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
# BoundaryConditions = "Neumann";
BoundaryConditions = "Periodic";

######
###### Choose one of the following levels of Laplacian discretization
######

# Order = "2";
# Order = "4";
# Order = "6";
Order = "8";

######
######
######


InitialConditions = VIni;
# InitialConditions = Sets.VIniTower;
# InitialConditions = Sets.CstStableTower(5.1, [0.5, 0.51]);
# InitialConditions = Sets.CstStableTowerRandom(5.0, [0.0, 1.0]);
# InitialConditions = Sets.CstStableMediumCstPerturb


ParameterSet = Set;
# ParameterSet = Sets.SetL2;
# ParameterSet = Sets.CstStable;