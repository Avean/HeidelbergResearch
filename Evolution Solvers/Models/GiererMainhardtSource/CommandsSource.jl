


# Perturbators 
ChangeStateLive(0.01,0.1,4.0,'S')
ChangeStateLive(0.55,0.6,4.0,'S')
ChangeStateLive(0.6,0.65,1.0,'S')
ChangeStateLive(0.7,0.75,1.0,'S')
ChangeStateLive(0.8,0.85,10.0,'S')

# Time step updates
UpdateDt(5e-1)
UpdateDt(1e-1)
UpdateDt(1e-2)
UpdateDt(1e-3)
UpdateDt(1e-4)
UpdateDt(1e-5)

# Set constant values
SetConstantLive(1.0,2.0, 1.0)
ChangeStateLive(0.01,0.1,1.001,'S')
ChangeStateLive(0.6,0.65,1.0,'S')


# Killers of heads
ChangeStateLive(0.3,0.8,0.01,'S')

# DDI One head + Wound
SetConstantLive(1.0,2.0,1.0)
ChangeStateLive(0.01,0.1,0.001,'P')
ChangeStateLive(0.6,0.65,0.7,'P')


# Wound
SetConstantLive(0.0,0.0,0.0)
ChangeStateLive(0.001,0.01,0.1,'P')
ChangeStateLive(0.6,0.65,0.1,'P')

# Different Diffusion
ChangeStateLive(0.01,0.05,5.1,'P')
ChangeStateLive(0.95,0.99,5.1,'P')
ChangeStateLive(0.65,0.7,5.1,'P')
ChangeStateLive(0.48,0.52,1.0,'P')

ChangeStateLive(0.65,0.7,5.1,'P')
ChangeStateLive(0.48,0.52,5.1,'P')