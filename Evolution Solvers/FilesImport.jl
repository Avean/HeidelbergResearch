
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
includet("BiffurcationPointsCompute.jl")

######
###### Viewer File
######

includet("Models/"*ModelName*"/"*ModelName*"Panel.jl")
includet("Models/"*ModelName*"/"*ModelName*"AsynViewer.jl")


