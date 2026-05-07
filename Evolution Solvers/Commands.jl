


# Perturbators 
ChangeStateLive(0.01,0.1,1.0,'S')
ChangeStateLive(0.55,0.6,4.0,'S')
ChangeStateLive(0.6,0.65,1.0,'S')
ChangeStateLive(0.7,0.75,1.0,'S')
ChangeStateLive(0.8,0.85,3.0,'S')

# Time step updates
UpdateDt(5e-1)
UpdateDt(1e-1)
UpdateDt(1e-2)
UpdateDt(1e-3)
UpdateDt(1e-4)
UpdateDt(1e-5)

# Set constant values
SetConstantLive(1.0,2.0)

# Killers of heads
ChangeStateLive(0.5,0.8,0.0,'S')