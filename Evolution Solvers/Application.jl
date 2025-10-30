using Revise
using LinearAlgebra
using LinearSolve
using SparseArrays
# using Plots

# using Distributed
# addprocs(1)  # dodaj jeden proces roboczy do obsługi viewer'a

using Base.Threads

# using PyCall

includet("FilesImport.jl")
include("Parametry.jl")

using .Solvers
using .Sets
using .Extractor
using .Viewer


# server_task =  Solvers.snapshot_server!()


SharedState.stop_simulation[] =true
sim_task = Threads.@spawn Solvers.run_simulation!(
    InitialConditions,
    ParameterSet,
    Scheme,
    BoundaryConditions,
    Order,
    dt,
    NonlinearityFunction
    )
    
    
    XVars = Viewer.setup_viewer(ParameterSet,dt)
    Viewer.viewer_loop!(ParameterSet, XVars)
    
########### NEW PART ###########
##



S2 = deepcopy(StructExtract(W));
display(norm(S2[1][:])/sqrt(SimParam.N))

# Hole = 250:600;
# S2[1][Hole] = 0.0.*ones(length(Hole));
Up = 490:510;
S2[1][Up] = 50.0.*ones(length(Up));

# S2[1][1:1000] += 15.0 .*(rand(SimParam.N) .- 1/2);

V = Iteration(VariablesVector(S2...),
            ParameterSet,
            100.0,
            Scheme,
            BoundaryCondsitions,
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
    plot!(subplot = i, Y, title = string(Fields[i]), ylims=(minimum(X) - 0.1, maximum(X) + 0.1));
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
