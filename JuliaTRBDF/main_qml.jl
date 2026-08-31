using Pkg

Pkg.activate(@__DIR__)

include(joinpath(@__DIR__, "src", "ReactionDiffusionApp.jl"))
include(joinpath(@__DIR__, "src", "QMLInterface.jl"))

using .ReactionDiffusionQML

ReactionDiffusionQML.run_qml_app(N = 300)
nothing


