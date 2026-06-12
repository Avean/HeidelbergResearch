module ρChange
    using ..SimParam
    using ..Nonlinearity
    using Observables

    ρObs = Observable(Nonlinearity.ρ)

    function UpdateRho(x1,x2, Val)
        x1 = round(Int, x1 / SimParam.dx)
        x2 = round(Int, x2 / SimParam.dx)
        Nonlinearity.ρ[x1:x2] .= Val; # Change production rate
        notify(ρObs)
    end


    function RhoSlope(Beg, End, ValB, ValE)
        i_start = max(1, round(Int, Beg * SimParam.N))
        i_end   = min(SimParam.N, round(Int, End * SimParam.N))

        Nonlinearity.ρ[i_start:i_end] .= range(
            ValB,
            ValE;
            length = i_end - i_start + 1
        )
        notify(ρObs)
    end 
end