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

    #######
    ####### Windows creation 
    #######
    Width = 1500
    F = Figure()
    WindowCreator.show_fig!(WindowCreator.screen2, F)

    Base.@kwdef mutable struct Bounds
        min::Float64 = 0.0
        max::Float64 = 1.0
    end

    # Slider
    Bν3 = Bounds()

    Sν3 =  Slider(F[2, 2], range = range(Bν3.min,Bν3.max,1000), startvalue = 0.5, width = Width)
    Lν3 = Label(F[1,2], @sprintf("Diffusion ν3 = %0.3f", Sν3.value[]))

    on(Sν3.value) do x
        # println("Slider moved to $(x)")
        Lν3.text[] = @sprintf("Diffusion ν3 = %0.3f", Sν3.value[])
    end

    BTν3Min = Textbox(F[2,1], placeholder = "Min")
    BTν3Max = Textbox(F[2,3], placeholder = "Max")

    BTν3Min.displayed_string[] = string(Bν3.min)
    BTν3Max.displayed_string[] = string(Bν3.max)

    BAν3 = Axis(F[3,2], width = Width, height = 20,     
                spinewidth = 0,           # usuwa ramkę
                xticksvisible = false,
                yticksvisible = false,
                xticklabelsvisible = false,
                yticklabelsvisible = false
                )



                
    on(BTν3Min.stored_string) do x
        try 
            Bν3.min = parse(Float64,x)

            xlims!(BAν3, Bν3.min, Bν3.max)
            Sν3.range = range(Bν3.min, Bν3.max, 1000)
                
            if Sν3.value[] < Bν3.min
                Sν3.value[] = Bν3.min
            end
        catch e
            display(Bν3.min)
            @warn "Wrong number in BTν3Min: $x" exception = e
        end
    end    

    on(BTν3Max.stored_string) do x
        try 
            Bν3.max = parse(Float64,x)

            xlims!(BAν3, Bν3.min, Bν3.max)
            Sν3.range = range(Bν3.min, Bν3.max, 1000)
                
            if Sν3.value[] > Bν3.max
                Sν3.value[] = Bν3.max
            end
        catch e
            display(Bν3.max)
            @warn "Wrong number in BTν3Max: $x" exception = e
        end
    end    


                
    BSν3 = scatter!(BAν3, [0.0, 0.1, 0.3, 0.5], [0.0, 0.0, 0.0, 0.0])
    xlims!(BAν3, Bν3.min, Bν3.max)
    ylims!(BAν3, -0.1, 0.1)

    
    for (x, y, label) in zip([0.0, 0.1, 0.3, 0.5], [0.0, 0.0, 0.0, 0.0].-0.05, ["A", "B", "C", "D"])
        text!(BAν3, label, position = (x, y), align = (:left, :center))
    end    
end