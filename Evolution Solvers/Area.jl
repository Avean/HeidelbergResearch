includet("FilesImport.jl")
include("Parametry.jl")

using .Nonlinearity
using Statistics

TP = range(0,ParameterSet.Coef.lbreak/ParameterSet.Coef.Slope,1000);
H1 = Nonlinearity.TimeSlope(ParameterSet, TP);
H2 = H1 .- 1 .- ParameterSet.Diff.D1 * (2pi)^2;
mean(H2)