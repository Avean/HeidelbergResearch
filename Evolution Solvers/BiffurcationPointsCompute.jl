module ComputeBifurctaionPoints

    using LinearAlgebra
    using ..Struktury
    using ..Sets
    

    function get_bif_points_linear(idx::Int; max_wn::Int = 2000)
        n = fieldcount(Struktury.VariablesVector)
        @assert 1 <= idx <= n "idx poza zakresem"

        Ds = Float64[]
        ks = Int[]

        D0 = [getfield(Sets.Par.Diff,i) for i in 1:fieldcount(Diffusions)]
        D0 = Diagonal(D0)
        D1 = copy(D0)

        for k in 1:max_wn
            wn2 = (k^2) * pi^2

            # D = 0 w wybranej współrzędnej
            
            D0[idx, idx] = 0.0

            # D = 1 w wybranej współrzędnej
            
            D1[idx, idx] = 1.0

            M0 = Sets.Jac .- wn2 .* D0
            M1 = Sets.Jac .- wn2 .* D1

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

        D0 = zeros(length(Ds))

        return Ds, D0, ks
    end

end