import sympy as sp
from scipy.optimize import fsolve
import numpy as np

from matrix_subsystems import get_principle_submatrices, get_stability_submatrices

gam = sp.Symbol('gam')
alp = sp.Symbol('alp')
alpp = sp.Symbol('alpp')
bet = sp.Symbol('bet')
lam = sp.Symbol('lam')

T = sp.Symbol('T')
B = sp.Symbol('B')
W = sp.Symbol('W')

c = 0.002
d = 0.01
k = 0.01 #0.01 #0.2
l = 0.01
p = 1.0 # p in [0,2]
q = 0.05
r = 0.35
s = 0.2 # 0.0
w = 0.001

alpha_ex = d / (k+w*p)
beta_ex = alpha_ex**2*c*p*s*q / (d**2*r)
gamma_ex = alpha_ex**3*c**2*p**2 / (d**3*r)
lambda_ex = alpha_ex*l / d
beta = [alpha_ex, beta_ex, gamma_ex, lambda_ex]

# Nonlinearities:
f1 = -alp*T + bet*(T*B)/(1+T*B) #alp*B - (1 - bet*B)*T
expr_T = (alp*B)/(1 - bet*B)

f2 = gam*B*B*W - (alpp + bet*T)*B

f3 = 1 - gam*B*B*W - lam*W
expr_W = 1/(gam*B*B + lam)

# General Jacobian matrix:
f1_T = sp.diff(f1, T)
f1_B = sp.diff(f1, B)
f1_W = sp.diff(f1, W)
f2_T = sp.diff(f2, T)
f2_B = sp.diff(f2, B)
f2_W = sp.diff(f2, W)
f3_T = sp.diff(f3, T)
f3_B = sp.diff(f3, B)
f3_W = sp.diff(f3, W)
Jac_gen = sp.Matrix([[f1_T, f1_B, f1_W],
                      [f2_T, f2_B, f2_W],
                      [f3_T, f3_B, f3_W]])

print("General Jacobian Matrix:")
sp.pprint(Jac_gen)
print("")

# Delta = sp.sqrt(gam*(gam - 4*(alp+bet)*alp*lam))
# # Delta = sp.Symbol('Delta', real=True)
# T2 = (-2*alp*bet*lam + gam + Delta) / (2*(bet**2*lam + gam))
# B2 = (gam + Delta) / (2*(alp + bet)*gam)
# W2 = (2*(alp+bet)*bet*lam + gam - Delta) / (2*lam*(bet**2*lam+gam))

# T3 = (-2*alp*bet*lam + gam - Delta) / (2*(bet**2*lam + gam))
# B3 = (gam - Delta) / (2*(alp + bet)*gam)
# W3 = (2*(alp+bet)*bet*lam + gam + Delta) / (2*lam*(bet**2*lam+gam))

# # det12 = f1_T*f2_B - f1_B*f2_T
# # det13 = f1_T*f3_W - f1_W*f3_T
# # det23 = f2_B*f3_W - f2_W*f3_B


# # det12_at_B2 = sp.simplify(det12.subs({T: T2, B: B2, W: W2}))
# # det13_at_B2 = sp.simplify(det13.subs({T: T2, B: B2, W: W2}))
# # det23_at_B2 = sp.simplify(det23.subs({T: T2, B: B2, W: W2}))

# # print("det12 at the equilibrium point 2:", det12_at_B2)
# # print("")
# # print("det13 at the equilibrium point 2:", det13_at_B2)
# # print("")
# # print("det23 at the equilibrium point 2:", det23_at_B2)
# # print("")
# # print("j11 at the equilibrium point 2:", sp.simplify(f1_T.subs({T: T2, B: B2, W: W2})))
# # print("")
# # print("j22 at the equilibrium point 2:", sp.simplify(f2_B.subs({T: T2, B: B2, W: W2})))
# # print("")

# # p = sp.Symbol('p')
# # sp.solve(sp.simplify(det12_at_B2.subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3]})) < 0, p)   # ((0.692478343306031 <= p) & (p < 2.31801672685591))
# # # sp.solve(sp.simplify(det13_at_B2.subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3]})) < 0, p)  # does not converge
# # sp.solve(sp.simplify(det23_at_B2.subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3]})) < 0, p)   # ((0.692478343306031 <= p) & (p < 0.702567944257956))

# # sp.solve(sp.simplify(f1_T.subs({T: T2, B: B2, W: W2, alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3]})) > 0, p)   # (-10.0 < p) & (p < -1.48936170212766)


# T2_val = T2.subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3]})
# B2_val = B2.subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3]})
# W2_val = W2.subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3]})
# T3_val = T3.subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3]})
# B3_val = B3.subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3]})
# W3_val = W3.subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3]})

# f2_eval = f2
# f2_eval = f2_eval.subs({W: expr_W, T: expr_T, alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3]})
# B0 = sp.solve(f2_eval, B)
# print(B0)
# B0_real = [np.real(np.complex128(b)) for b in B0 if abs(sp.im(b)) < 10e-15]

# if B0_real:
#     if max(B0_real) > 1.0e5:
#         # Take second largest value if the largest is unreasonably high
#         B0_real_sorted = sorted(B0_real)
#         B0_max = B0_real_sorted[-2]
#     else:
#         B0_max = max(B0_real)
# else:
#     print("Could not find real steady state!")
#     B0_max = 42 # Some random value


# T0 = expr_T.subs({B: B0_max, alp: beta[0], bet: beta[1]})
# W0 = expr_W.subs({B: B0_max, gam: beta[2], lam: beta[3]})

# # test steady state:
# print("Testing steady state:")
# print(f"Steady state values: T={T0}, B={B0_max}, W={W0}")
# print(f1.subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3], T: T0, B: B0_max, W: W0}))
# print(f2.subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3], T: T0, B: B0_max, W: W0}))
# print(f3.subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3], T: T0, B: B0_max, W: W0}))
# print("")

# print("Compare to paper steady state values:")
# print(f"Steady state values: T={T2_val}, B={B2_val}, W={W2_val}")
# print(f"Steady state values: T={T3_val}, B={B3_val}, W={W3_val}")

# j11 = sp.diff(f1, T).subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3], T: T0, B: B0_max, W: W0})
# j12 = sp.diff(f1, B).subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3], T: T0, B: B0_max, W: W0})
# j13 = sp.diff(f1, W).subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3], T: T0, B: B0_max, W: W0})
# j21 = sp.diff(f2, T).subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3], T: T0, B: B0_max, W: W0})
# j22 = sp.diff(f2, B).subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3], T: T0, B: B0_max, W: W0})
# j23 = sp.diff(f2, W).subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3], T: T0, B: B0_max, W: W0})
# j31 = sp.diff(f3, T).subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3], T: T0, B: B0_max, W: W0})
# j32 = sp.diff(f3, B).subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3], T: T0, B: B0_max, W: W0})
# j33 = sp.diff(f3, W).subs({alp: beta[0], bet: beta[1], gam: beta[2], lam: beta[3], T: T0, B: B0_max, W: W0})

# Jacobian = sp.Matrix([[j11, j12, j13],
#                       [j21, j22, j23],
#                       [j31, j32, j33]])
# print("Jacobian Matrix at the equilibrium point:")
# sp.pprint(Jacobian)
# print("")
# print("Eigenvalues of the Jacobian Matrix:")
# eigenvals = Jacobian.eigenvals()
# for val, mult in eigenvals.items():
#     print(f"Eigenvalue: {val}, Multiplicity: {mult}")
# print("")
# print("Jacobian Matrix at the equilibrium point for copying:")
# print(Jacobian)
# Jac_np = np.array(Jacobian).astype(np.float64)
# print('Check Jacobian for unstable subsystems:')
# principal_submatrices = get_principle_submatrices(Jac_np)
# submatrices_stable = []
# for dim, submats in principal_submatrices.items():
#     print('')
#     print(f"Dimension: {dim}x{dim}")
#     print('Stability:')
#     submatrices_stable.append(get_stability_submatrices(submats))
#     # print('All submatrices:')
#     # for submatrix in submats:
#     #     print(submatrix)
# print('')
# if submatrices_stable == [True for dim in principal_submatrices.keys()]:
#     print('All principal submatrices are stable!')
# else:
#     print('There is at least one unstable principal submatrix.')

# # Jac_variant_2 = np.array([[-1.54835418090850, -0.968673468738539, 0, -0.519279522931973, 1.00000000000000], [-0.0713698848523722, -1, 0, 0, 0], [1.15960000000000, 0, -1, 0, 11.5964000000000], [-0.896792686154216, 0, 0.0760938750398222, -1, 0], [1, 0, 0, 0, -1]])
