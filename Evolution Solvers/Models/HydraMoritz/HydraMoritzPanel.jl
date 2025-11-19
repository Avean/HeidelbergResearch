module WindowCreator
   
    using GLMakie

    # Screen References -> 2 windows, initialy no windows
    screen1 = Ref{Union{Nothing, GLMakie.Screen}}(nothing)
    screen2 = Ref{Union{Nothing, GLMakie.Screen}}(nothing)

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

    export ResetPanel

    #######
    ####### Windows creation 
    #######
    Width = 1500
    F = Figure()
    WindowCreator.show_fig!(WindowCreator.screen2, F)

    Base.@kwdef mutable struct Bounds
        min::Float64
        max::Float64
    end

    BmM = Vector(undef,fieldcount(Diffusions))    # Bound min max
    SL = Vector(undef,fieldcount(Diffusions))     # Slider
    LTit = Vector(undef,fieldcount(Diffusions))   # Title with caption
    TBmin = Vector(undef,fieldcount(Diffusions))  # Textbox min
    TBmax = Vector(undef,fieldcount(Diffusions))  # Textbox max
    AX = Vector(undef,fieldcount(Diffusions))     # Axis Vector

    for i in 1:fieldcount(Diffusions)

        BmM[i] = Bounds(0.0, 2*getfield(Sets.Par.Diff,i)+1e-6)

        SL[i] =  Slider(F[2+4*(i-1), 2], range = range(BmM[i].min, BmM[i].max, 1000), startvalue = getfield(Sets.Par.Diff,i), width = Width)
        LTit[i] = Label(F[1+4*(i-1),2], @sprintf("Diffusion %s = %0.3g", string(fieldnames(Diffusions)[i]),SL[i].value[]))

        on(SL[i].value) do x
            LTit[i].text[] = @sprintf("Diffusion %s = %0.3g", string(fieldnames(Diffusions)[i]), SL[i].value[])
            setfield!(Sets.Par.Diff,i,SL[i].value[])
        end

        TBmin[i] = Textbox(F[2+4*(i-1),1], placeholder = "Min")
        TBmax[i] = Textbox(F[2+4*(i-1),3], placeholder = "Max")

        TBmin[i].displayed_string[] = string(BmM[i].min)
        TBmax[i].displayed_string[] = string(BmM[i].max)

        AX[i] = Axis(F[3+4*(i-1),2], width = Width, height = 20,     
                    spinewidth = 0,           # usuwa ramkę
                    xticksvisible = false,
                    yticksvisible = false,
                    xticklabelsvisible = false,
                    yticklabelsvisible = false
                    )
                    
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
        end    

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
        end    
                    
        xlims!(AX[i], BmM[i].min, BmM[i].max)
        
        
        for (x, y, label) in zip(Sets.BifurcationData.PointsX[i], Sets.BifurcationData.Points0[i].-0.05, Sets.BifurcationData.Names[i])
            text!(AX[i], label, position = (x, y), align = (:left, :center))
        end
        scatter!(AX[i], Sets.BifurcationData.PointsX[i],Sets.BifurcationData.Points0[i])    
        ylims!(AX[i], -0.1, 0.1)
    end


    function ResetPanel()
        Sets.ResetParameters()
        
        for i in 1:fieldcount(Diffusions)
                   
            BmM[i] = Bounds(0.0, 2*getfield(Sets.Par.Diff,i)+1e-6)
            SL[i].range[] = range(BmM[i].min, BmM[i].max, 1000)
            
            Sets.ResetParameters()
            
            TBmin[i].displayed_string[] = string(BmM[i].min)
            TBmax[i].displayed_string[] = string(BmM[i].max)
            
            set_close_to!(SL[i], getfield(Sets.Par.Diff,i))
            
            xlims!(AX[i], BmM[i].min, BmM[i].max)
        end
    end
end

