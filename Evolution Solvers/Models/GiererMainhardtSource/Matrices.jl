using LinearAlgebra
using Random
using Plots

Random.seed!(1234)

# Losowa macierz 3x3, ale z trzema dodatnimi wartościami własnymi
sigma = 5.0

λ0 = [-0.5, 1.0, 2.0]
Λ = Diagonal(λ0)

P = randn(3, 3)

while abs(det(P)) < 1e-6
    P = randn(3, 3)
end

A = P * Λ * inv(P)

println("Macierz A:")
println(A)

println("Początkowe wartości własne macierzy A:")
println(eigvals(A))

# Który wiersz skalujemy
row_to_scale = 2

# Dodatnie wartości τ
taus = range(0.0, 5.0, length = 600)

tol = 1e-8

τ_real = Float64[]
λ_real = Float64[]

τ_complex = Float64[]
λ_complex = Float64[]

for τ in taus
    B = copy(A)

    # Skalowanie wybranego wiersza przez τ
    B[row_to_scale, :] .*= τ

    eigs = eigvals(B)

    for λ in eigs
        if abs(imag(λ)) < tol
            push!(τ_real, τ)
            push!(λ_real, real(λ))
        else
            push!(τ_complex, τ)
            push!(λ_complex, real(λ))
        end
    end
end

scatter(
    τ_real,
    λ_real,
    color = :blue,
    label = "wartość własna rzeczywista",
    markersize = 3,
    markerstrokewidth = 0,
    xlabel = "τ",
    ylabel = "Re(λ)",
    title = "Części rzeczywiste wartości własnych po skalowaniu wiersza"
)

scatter!(
    τ_complex,
    λ_complex,
    color = :green,
    label = "wartość własna zespolona",
    markersize = 3,
    markerstrokewidth = 0
)

display(current())