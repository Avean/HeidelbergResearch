using LinearAlgebra
using Printf

"""
    get_bif_points_linear(J, Ddiag, idx; max_wn=200)

Zakładamy, że dla ustalonego k wyznacznik
f(D) = det(J - k^2 * π^2 * D(D))
jest liniowy względem D = Ddiag[idx].

J     - macierz Jacobiego (n×n)
Ddiag - wektor diagonalny D (długości n)
idx   - indeks (1-based) elementu diagonali traktowanego jako zmienna D
max_wn - maksymalne k do sprawdzenia

Zwraca: (Ds, ks)
- Ds: dodatnie wartości D(k), dla których det = 0
- ks: odpowiadające im k
"""
function get_bif_points_linear(J::AbstractMatrix, Ddiag::AbstractVector, idx::Int; max_wn::Int = 200)
    n = length(Ddiag)
    @assert size(J,1) == size(J,2) == n "J musi być macierzą kwadratową zgodną z Ddiag"
    @assert 1 <= idx <= n "idx poza zakresem"

    Ds = Float64[]
    ks = Int[]

    for k in 1:max_wn
        wn2 = (k^2) * pi^2

        # D = 0 w wybranej współrzędnej
        D0 = copy(Ddiag)
        D0[idx] = 0.0

        # D = 1 w wybranej współrzędnej
        D1 = copy(Ddiag)
        D1[idx] = 1.0

        M0 = J .- wn2 .* Diagonal(D0)
        M1 = J .- wn2 .* Diagonal(D1)

        f0 = det(M0)  # f(0) = b
        f1 = det(M1)  # f(1) = a + b

        a = f1 - f0
        b = f0

        # jeśli prawie brak zależności od D – pomijamy
        if abs(a) < 1e-12
            continue
        end

        Dstar = -b / a

        if Dstar > 0
            push!(Ds, Dstar)
            push!(ks, k)
        else
            # jak chcesz, możesz zamiast break dać continue,
            # ale w duchu Twojego poprzedniego kodu wychodzimy z pętli
            break
        end
    end

    return Ds, ks
end

Jac = [
    -1.0865168027146117   -0.05283439313362805   0.0                 -0.028643941383072864   0.9999999999999997
    -3.2510268461888066   -1.0                   0.0                  0.0                    0.0
    44.1374534855984       0.0                  -1.0                  0.0                   44.1374534855984
    -11.026231669288988    0.0                   0.5135733514383242  -1.0                    0.0
    1.0                    0.0                   0.0                  0.0                   -1.0
]




function print_bif_table(Ds, ks)
    println("===========================================")
    println("       PUNKTY BIFURKACJI (D*, k)")
    println("===========================================")
    for (i, (k, Dval)) in enumerate(zip(ks, Ds))
        println("Punkt $(i):")
        @printf("   k       = %d\n", k)
        @printf("   Wartość = %.12g\n", Dval)
        println("-------------------------------------------")
    end
end

Ddiag = [0.0, 3.8154e-5, 0.4433, 6.0713e-3, 4e-4]
idx = 3   # modyfikujemy 3. współrzędną (jak Twoje d w Pythonie)
Ds, ks = @time get_bif_points_linear(Jac, Ddiag, idx; max_wn=20000)
print_bif_table(Ds, ks)







