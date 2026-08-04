#!/usr/bin/env julia

using Pkg

# External packages imported by src/ReactionDiffusionApp.jl.
# Julia standard-library modules such as SparseArrays, LinearAlgebra,
# Random, Printf, Statistics, and Pkg do not need to be installed.
const REQUIRED_PACKAGES = [
    "CxxWrap",
    "CairoMakie",
    "GLMakie",
    "DifferentialEquations",
    "OrdinaryDiffEqSDIRK",
    "LaTeXStrings",
    "QML",
    "QMLMakie",
]


function install_dependencies()
    project_dir = @__DIR__

    println("Julia version: ", VERSION)
    println("Activating project: ", project_dir)

    Pkg.activate(project_dir)

    declared_packages = Pkg.project().dependencies
    missing_packages = filter(
        package_name -> !haskey(declared_packages, package_name),
        REQUIRED_PACKAGES,
    )

    if isempty(missing_packages)
        println("All required external packages are already declared.")
    else
        println(
            "Adding missing packages: ",
            join(missing_packages, ", "),
        )

        Pkg.add(missing_packages)
    end

    println("Instantiating the project environment...")
    Pkg.instantiate()

    println("Precompiling dependencies...")
    Pkg.precompile()

    println()
    println("Dependency installation completed successfully.")
    println("Start the application with:")
    println("    julia --threads=auto --project=. main_qml.jl")
    println()
    println("The previous GLMakie interface remains available with:")
    println("    julia --threads=auto --project=. main.jl")

    return nothing
end


if abspath(PROGRAM_FILE) == @__FILE__
    install_dependencies()
end
