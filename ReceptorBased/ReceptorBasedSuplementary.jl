using ..Sets

eta1 = Sets.SetL1.Coef.m1 / Sets.SetL1.Coef.u1;
eta2 = Sets.SetL1.Coef.m2 / Sets.SetL1.Coef.u2;
eta3 = Sets.SetL1.Coef.m3 / Sets.SetL1.Coef.u3;
m1 = Sets.SetL1.Coef.m1;
m2 = Sets.SetL1.Coef.m2;
m3 = Sets.SetL1.Coef.m3;

mul = Sets.SetL1.Coef.u2; 

v_p = (m2*eta1 + eta3)/(2*eta1*(mul + eta3)) + sqrt((m2*eta1 - eta3)^2 - 4*m2*mul*eta1)/(2*eta1*(mul + eta3))
u_p = eta1 - 1/v_p
w_p = eta3 - eta3/eta1/v_p

v_m = (m2*eta1 + eta3)/(2*eta1*(mul + eta3)) - sqrt((m2*eta1 - eta3)^2 - 4*m2*mul*eta1)/(2*eta1*(mul + eta3))
u_m = eta1 - 1/v_m
w_m = eta3 - eta3/eta1/v_m