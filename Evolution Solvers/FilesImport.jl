# Choose a model by selecting a folder

ModelName = "HydraMoritz";
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

######
###### Viewer File
######

includet("Models/"*ModelName*"/"*ModelName*"AsynViewer.jl")



