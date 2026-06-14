# main.jl

using Pkg

# Activate the project located in the same folder as this file.
Pkg.activate(@__DIR__)

# Load the application module.
include(joinpath(@__DIR__, "src", "ReactionDiffusionApp.jl"))

using .ReactionDiffusionApp

# Start the application.
ReactionDiffusionApp.run_app();

1