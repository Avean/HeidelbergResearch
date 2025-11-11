# Choose a model by selecting a folder

# ModelName = "HydraDietmar";
# ModelName = "Test1";
ModelName = "GiereMeinhardt";
# ModelName = "ReceptorBased";

######
###### Include the model files
######

includet("Models/"*ModelName*"/"*ModelName*"Modules.jl")

######
###### Solvers File
######

includet("DiffusionMatrices.jl")
includet("SolverSteps.jl")
includet("AsynViewer.jl")


