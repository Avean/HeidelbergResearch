using Revise
using LinearAlgebra
using LinearSolve
using SparseArrays
# using Plots

using Base.Threads

# Choose a model by selecting a folder

# ModelName = "HydraDietmarA";  # Linear WntDiff with WntLoc coupling only 
ModelName = "HydraDietmarB";  # Linear WntDiff with WntLoc coupling only 
# ModelName = "HydraMoritz5eq";   # Original model with 5 equations   
# ModelName = "HydraMoritz5eqA";  # Linear WntDiff with SD and WntLoc coupling
# ModelName = "HydraMoritz5eqB";  # Linear WntDiff with WntLoc coupling only 
# ModelName = "HydraMoritz4eqA";  # 4 equations without SD, original parameters
# ModelName = "HydraMoritz4eqB";  # 4 equations without SD, changed paramters
# ModelName = "Test1";
# ModelName = "GiereMeinhardt";
# ModelName = "ReceptorBased";



includet("FilesImport.jl")

using .Solvers
using .Sets: DisplayDDI
using .Extractor
using .Viewer
using .Panel: SetPanel, ResetPanel

# server_task =  Solvers.snapshot_server!()

@time SetPanel()

XVars = Viewer.setup_viewer();



# include("Parametry.jl")

Viewer.stop_simulation!(XVars)
Viewer_task = @async begin
    Viewer.viewer_loop!(XVars)
end
sim_task = Threads.@spawn Solvers.run_simulation!(
    Sets.Ini,
    Scheme,
    BoundaryConditions,
    Order,
    NonlinearityFunction
    )




using GLMakie
using .SharedState
F = Figure()
ax = Axis(F[1, 1])
Y = frame_buffer[][1];
lines!(ax, Y.WntDiff)
# lines!(ax, Y.DkkC)

@time begin
    Y = UnpackStruct(Sets.Ini);
    A = hcat(Y...);
    B = eachrow(A);
    Sets.Fun.(B, Ref(Sets.β));
end
    
# @async Viewer.RecordAnimation(20.0, "HydraStable1.9.mp4", 0.0)


SharedState.pause_simulation[] = true
SharedState.pause_simulation[] = false


##


########### NEW PART ###########
##



# S2 = deepcopy(StructExtract(W));
# display(norm(S2[1][:])/sqrt(SimParam.N))

# Hole = 250:600;
# S2[1][Hole] = 0.0.*ones(length(Hole));
# Up = 490:510;
# S2[1][Up] = 50.0.*ones(length(Up));

# # S2[1][1:1000] += 15.0 .*(rand(SimParam.N) .- 1/2);

# V = Iteration(VariablesVector(S2...),
#             ParameterSet,
#             100.0,
#             Scheme,
#             BoundaryCondsitions,
#             Order,
#             dt,
#             NonlinearityFunction);

# print("Done")

# Fields = fieldnames(VariablesVector);
# FieldsNum = length(Fields);
# NFields = 1:FieldsNum;

# P = plot(layout = (FieldsNum,1));
# for i in NFields
#     X = getfield(V, Fields[i]);
#     Y = getfield(W, Fields[i]);
#     plot!(subplot = i, X, title = string(Fields[i]), ylims=(minimum(X) - 0.1, maximum(X) + 0.1));
#     plot!(subplot = i, Y, title = string(Fields[i]), ylims=(minimum(X) - 0.1, maximum(X) + 0.1));
# end
# display(P)

##
########### NEW PART ###########

# Um = zeros(1000);

# Vm = 0.0.*ones(1000);
# Vm[1:100] = ones(100);
# # Vm[901:1000] = ones(100);

# Wm = 0.0.*ones(1000);
# Wm[1:100] = ones(100);
# Wm[901:1000] = ones(100);
