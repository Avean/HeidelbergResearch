using Revise

includet("FilesImport.jl")
include("Parametry.jl")

using .Viewer


XVars = setup_viewer(ParameterSet,dt)
Viewer.server_loop!(ParameterSet, XVars)