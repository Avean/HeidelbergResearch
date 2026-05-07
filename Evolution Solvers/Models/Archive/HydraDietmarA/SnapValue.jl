using JLD2

DV = Vector{Float64}()
κV = Vector{Float64}()
UV = Vector{Vector{Float64}}()




κ = Sets.Par.Coef.κ
D = Sets.Par.Diff.D1
Un = SharedState.frame_buffer[][1].u

push!(DV, D)
push!(κV, κ)
push!(UV, Un)



@save "SnapValue.jld2" DV κV UV