using ..Struktury
using ..SimParam
using .Sets


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

######
###### Choosing Nonlinearity
######

# NonlinearityFunction = "Nonlinearity 1"; # Classical model
NonlinearityFunction = "Nonlinearity 2"; # Reduced 3eq model
######
###### Initial conditions and parameters
######


