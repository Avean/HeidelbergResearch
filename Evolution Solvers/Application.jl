using Revise
using LinearAlgebra
using LinearSolve
using SparseArrays
using Plots


# using PyCall

# Choose a model by selecting a folder

ModelName = "HydraDietmar";
# ModelName = "Test1";
# ModelName = "GiereMeinhardt";
# ModelName = "ReceptorBased";

######
###### Include the model files
######

includet("Models/"*ModelName*"/"*ModelName*"Modules.jl")
include("Models/"*ModelName*"/"*ModelName*"Variables.jl")

######
###### Solvers File
######

includet("DiffusionMatrices.jl")
includet("SolverSteps.jl")

using .Solvers
using .Sets
using .Extractor

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


# InitialConditions = VIni;
# InitialConditions = Sets.VIniTower;
InitialConditions = Sets.CstStableTower(5.1, [0.5, 0.51]);
# InitialConditions = Sets.CstStableTowerRandom(5.0, [0.0, 1.0]);
# InitialConditions = Sets.CstStableMediumCstPerturb


ParameterSet = Set;
# ParameterSet = Sets.SetL2;
# ParameterSet = Sets.CstStable;

W = Iteration(InitialConditions,
            ParameterSet,
            10.0,
            Scheme,
            BoundaryConditions,
            Order,
            dt,
            NonlinearityFunction);

########### NEW PART ###########
##



S2 = deepcopy(StructExtract(W));
display(norm(S2[1][:])/sqrt(SimParam.N))

# Hole = 250:600;
# S2[1][Hole] = 0.0.*ones(length(Hole));
Up = 500:510;
S2[1][Up] = 50.0.*ones(length(Up));

# S2[1][1:1000] += 15.0 .*(rand(SimParam.N) .- 1/2);

V = Iteration(VariablesVector(S2...),
            ParameterSet,
            100.0,
            Scheme,
            BoundaryConditions,
            Order,
            dt,
            NonlinearityFunction);

print("Done")

Fields = fieldnames(VariablesVector);
FieldsNum = length(Fields);
NFields = 1:FieldsNum;

P = plot(layout = (FieldsNum,1));
for i in NFields
    X = getfield(V, Fields[i]);
    Y = getfield(W, Fields[i]);
    plot!(subplot = i, X, title = string(Fields[i]), ylims=(minimum(X) - 0.1, maximum(X) + 0.1));
    plot!(subplot = i, X, title = string(Fields[i]), ylims=(minimum(X) - 0.1, maximum(X) + 0.1));
end
display(P)

##
########### NEW PART ###########

# Um = zeros(1000);

# Vm = 0.0.*ones(1000);
# Vm[1:100] = ones(100);
# # Vm[901:1000] = ones(100);

# Wm = 0.0.*ones(1000);
# Wm[1:100] = ones(100);
# Wm[901:1000] = ones(100);
