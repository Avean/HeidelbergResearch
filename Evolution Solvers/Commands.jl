


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
SetConstantLive(1.0,2.0)
ChangeStateLive(0.01,0.1,1.001,'S')
ChangeStateLive(0.6,0.65,1.0,'S')


# Killers of heads
ChangeStateLive(0.001,0.2,0.00,'S')

# DDI One head + Wound
SetConstantLive(1.0,2.0)
ChangeStateLive(0.01,0.1,0.01,'P')
ChangeStateLive(0.6,0.65,0.7,'P')


# Wound
SetConstantLive(0.0,0.0)
ChangeStateLive(0.001,0.01,0.1,'P')
ChangeStateLive(0.6,0.65,0.1,'P')
ChangeStateLive(0.7,0.75,0.5,'P')


# Different Diffusion
ChangeStateLive(0.01,0.05,5.1,'P')
ChangeStateLive(0.95,0.99,5.1,'P')
ChangeStateLive(0.65,0.7,1.1,'P')
ChangeStateLive(0.48,0.52,1.0,'P')

ChangeStateLive(0.65,0.7,5.1,'P')
ChangeStateLive(0.48,0.52,5.1,'P')


Nonlinearity.UpdateRho(0.1,0.9,1.0)
# Play with different ρ

using ..ρChange
ρChange.RhoSlope(0.0,1.0,1.0,1.0)

# Scenarion 1. Alignment of parts
ρChange.RhoSlope(0.0,1.0,0.8,1.2)

SetConstantLive(1.0,2.0)

SetConstantLive(0.0,0.0)

ChangeStateLive(0.001,0.01,0.5,'P')
ChangeStateLive(0.99,0.999,0.5,'P')
ChangeStateLive(0.50,0.51,0.5,'P')


UpdateDt(1e-5)
UpdateDt(1e-2)
UpdateDt(5e-1)


# Scenario 2, Shifting parts
ρChange.RhoSlope(0.0,0.5,1.0,1.3)
ρChange.RhoSlope(0.5,1.0,0.7,1.0)

SetConstantLive(1.0,2.0)

SetConstantLive(0.0,0.0)

ChangeStateLive(0.001,0.01,0.5,'P')
ChangeStateLive(0.99,0.999,0.5,'P')
ChangeStateLive(0.50,0.51,0.5,'P')



UpdateDt(1e-5)
UpdateDt(1e-2)
UpdateDt(5e-1)


# Scenario 3, Anti Alignment foot Glue
ρChange.RhoSlope(0.0,0.5,1.0,0.7)
ρChange.RhoSlope(0.5,1.0,1.0,1.3)

SetConstantLive(1.0,2.0)

SetConstantLive(0.0,0.0)

ChangeStateLive(0.001,0.01,0.5,'P')
ChangeStateLive(0.99,0.999,0.5,'P')
ChangeStateLive(0.50,0.51,0.5,'P')



UpdateDt(1e-5)
UpdateDt(1e-2)
UpdateDt(5e-1)