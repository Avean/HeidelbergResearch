module WindowCreator
   
    using GLMakie

    # Screen References -> 2 windows, initialy no windows
    screen1 = Ref{Union{Nothing, GLMakie.Screen}}(nothing)
    screen2 = Ref{Union{Nothing, GLMakie.Screen}}(nothing)
    screen3 = Ref{Union{Nothing, GLMakie.Screen}}(nothing)

    function show_fig!(screen_ref::Base.RefValue, fig::Figure) # Creating new window
        if screen_ref[] === nothing || !isopen(screen_ref[])
            # New Screen if no screen or closed
            screen_ref[] = display(GLMakie.Screen(), fig)
        else #Update existing screen
            screen_ref[] = display(screen_ref[], fig)
            return
        end
    end
end


module Panel

    using GLMakie
    using Printf
    using ..WindowCreator
    using ..Struktury
    using ..Sets
    using ..Solvers
    using ..ComputeBifurctaionPoints: get_bif_points_linear, CheckDDI

    

    #######
    ####### Windows creation 
    #######
    Width = 1500


    Base.@kwdef mutable struct Bounds
        min::Float64
        max::Float64
    end

    BmM = Vector(undef,fieldcount(Diffusions))    # Bound min max
    SL = Vector(undef,fieldcount(Diffusions))     # Slider
    LTit = Vector(undef,fieldcount(Diffusions))   # Title with caption
    LBif = Vector(undef,fieldcount(Diffusions))   # Title with bifurcation
    TBmin = Vector(undef,fieldcount(Diffusions))  # Textbox min
    TBmax = Vector(undef,fieldcount(Diffusions))  # Textbox max
    AX = Vector(undef,fieldcount(Diffusions))     # Axis Vector

    ODs = [Observable(Float64[]) for i in 1:fieldcount(Diffusions)] 
    OKs = [Observable(String[]) for i in 1:fieldcount(Diffusions)]
    ODDI = Observable("Initializing...")

    function SetPanel()
        F = Figure()
        WindowCreator.show_fig!(WindowCreator.screen2, F)
        
        for i in 1:fieldcount(Diffusions)

            BmM[i] = Bounds(0.0, max(2*getfield(Sets.Par.Diff,i),1e-8))

            SL[i] =  Slider(F[3+5*(i-1), 2], range = range(BmM[i].min, BmM[i].max, 1000), startvalue = getfield(Sets.Par.Diff,i), width = Width)
            LTit[i] = Label(F[1+5*(i-1),2], @sprintf("Diffusion %s = %0.3g", string(fieldnames(VariablesVector)[i]),SL[i].value[]))


            ############# Slider On Change #############
            on(SL[i].value) do x
                LTit[i].text[] = @sprintf("Diffusion %s = %0.3g", string(fieldnames(VariablesVector)[i]), SL[i].value[])
                setfield!(Sets.Par.Diff,i,SL[i].value[])
                Solvers.UpdateDiffMatrix()
                for i in 1:fieldcount(Diffusions)
                    ODs[i][], OKs[i][] =UpdateBiffPlot(i, AX[i], BmM[i], LBif[i])
                    notify(ODs[i])
                    notify(OKs[i])
                end
                ODDI[] = CheckDDI()
            end
            ############# Slider On Change #############
            
            TBmin[i] = Textbox(F[3+5*(i-1),1], placeholder = "Min")
            TBmax[i] = Textbox(F[3+5*(i-1),3], placeholder = "Max")
            
            TBmin[i].displayed_string[] = @sprintf("%0.3g",BmM[i].min)
            TBmax[i].displayed_string[] = @sprintf("%0.3g",BmM[i].max)
            
            AX[i] = Axis(F[4+5*(i-1),2], width = Width, height = 20,     
            spinewidth = 0,           # usuwa ramkę
            xticksvisible = false,
            yticksvisible = false,
            xticklabelsvisible = false,
            yticklabelsvisible = false
            )
            
            ####### on Texbox change #######
            ########### Min ################
            on(TBmin[i].stored_string) do x
                try 
                    BmM[i].min = parse(Float64,x)
                    
                    xlims!(AX[i], BmM[i].min, BmM[i].max)
                    SL[i].range[] = range(BmM[i].min, BmM[i].max, 1000)
                    
                    if SL[i].value[] < BmM[i].min
                        set_close_to!(SL[i], BmM[i].min)
                    end
                catch e
                    display(BmM[i].min)
                    @warn "Wrong number in TBmin: $x" exception = e
                end

                for i in 1:fieldcount(Diffusions)
                    UpdateBiffPlot(i, AX[i], BmM[i], LBif[i])
                end
                ODDI[] = CheckDDI()
            end    

            ########### Max ################
            on(TBmax[i].stored_string) do x
                try 
                    BmM[i].max = parse(Float64,x)

                    xlims!(AX[i], BmM[i].min, BmM[i].max)
                    SL[i].range[] = range(BmM[i].min, BmM[i].max, 1000)
                        
                    if SL[i].value[] > BmM[i].max
                        set_close_to!(SL[i], BmM[i].max)
                    end
                catch e
                    display(BmM[i].max)
                    @warn "Wrong number in TBmax: $x" exception = e
                end

                for i in 1:fieldcount(Diffusions)
                    UpdateBiffPlot(i, AX[i], BmM[i], LBif[i])
                end
                ODDI[] = CheckDDI()
            end    
            ####### on Texbox change #######
                        


            xlims!(AX[i], BmM[i].min, BmM[i].max)

            LBif[i] = Label(F[2+5*(i-1),2], " ")
            
            # UpdateBiffPlot(i, AX[i], BmM[i], LBif[i])
            
            ODs0 = ODs[i]

            Pos = @lift(Point2f.($ODs0, -0.05))
            
            text!(AX[i], OKs[i], position = Pos, align = (:left, :center))
            scatter!(AX[i], ODs[i], @lift(0.0 .* $ODs0))
            
            ODs[i][], OKs[i][] = UpdateBiffPlot(i, AX[i], BmM[i], LBif[i])   

            ylims!(AX[i], -0.1, 0.1)
        end
        ODDI[] = CheckDDI()

        display(ODDI)
        Label(F[5*fieldcount(Diffusions),2], ODDI)

    end



    function UpdateBiffPlot(i::Int64, AX::Axis, BM, LBif)
        Ds, Ks = get_bif_points_linear(i)
   
        Kt = "k".*string.(Ks)
        
        if !isempty(Ds)
            L = length(Ds)
            DM = [minimum(Ds), maximum(Ds)]
            DK = [Kt[findmin(Ds)[2]], Kt[findmax(Ds)[2]]]

            LBif.text[] = @sprintf("Total number of bifurcation points = %d               [Min %s for %0.3g       Max  %s for %.3g]", L, string(DK[1]), DM[1],  string(DK[2]),DM[2])
        else
            LBif.text[] = "No bifurcation points found"
        end


        # return FilterBifSet(BM, Ds, Kt)
        return (Ds, Kt)
    end

    function FilterBifSet(BM, Ds, Kt)
        
        NumberOfBins = 100
        CapacityPerBin = 10

        DsBins = [Float64[] for i in 1:NumberOfBins]
        KtBins = [String[] for i in 1:NumberOfBins]
        
        dM = (BM.max - BM.min)/NumberOfBins
        
        for i in eachindex(Ds)
            if Ds[i] >= BM.min && Ds[i] <= BM.max
                bin_index = Int(floor((Ds[i] - BM.min)/dM)) + 1

                if length(DsBins[bin_index]) < CapacityPerBin
                    push!(DsBins[bin_index], Ds[i])
                    push!(KtBins[bin_index], Kt[i])
                end
            end
        end
        
        Dsf = reduce(vcat, DsBins)
        Ktf = reduce(vcat, KtBins)

        return Dsf, Ktf
    end

    function ResetPanel()
        Sets.ResetParameters()
        
        for i in 1:fieldcount(Diffusions)
                   
            BmM[i] = Bounds(0.0, max(2*getfield(Sets.Par.Diff,i),1e-8))
            SL[i].range[] = range(BmM[i].min, BmM[i].max, 1000)
            
            Sets.ResetParameters()
            
            TBmin[i].displayed_string[] = string(BmM[i].min)
            TBmax[i].displayed_string[] = string(BmM[i].max)
            
            set_close_to!(SL[i], getfield(Sets.Par.Diff,i))
            
            xlims!(AX[i], BmM[i].min, BmM[i].max)
        end
    end
end

