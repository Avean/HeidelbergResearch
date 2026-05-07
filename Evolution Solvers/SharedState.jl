module SharedState
    export request_frame, frame_buffer, modify_state, ModifyState

    const request_frame = Ref(false)
    const frame_buffer  = Ref{Any}(nothing)
    const stop_simulation = Ref(false)
    const pause_simulation = Ref(false)
    const stop_viewer = Ref(false)
    const modify_state = Ref(false)
    const modify_state_type = Ref{Char}('P') # 'S' for Set / 'P' for Perturb

    function ModifyState(Type::Char) # Function to modify curent state (S - Set / P - Perturb)
        modify_state[] = true
        modify_state_type[] = Type
    end

end