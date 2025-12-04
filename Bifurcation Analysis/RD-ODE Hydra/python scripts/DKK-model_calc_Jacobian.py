import sympy as sp
from scipy.optimize import fsolve
import numpy as np

from matrix_subsystems import get_principle_submatrices, get_stability_submatrices

Wl = sp.Symbol('Wl')
A = sp.Symbol('A')
Wd = sp.Symbol('Wd')
C = sp.Symbol('C')
S = sp.Symbol('S')
b1 = sp.Symbol('b1')
b2 = sp.Symbol('b2')
b3 = sp.Symbol('b3')
b4 = sp.Symbol('b4')
b5 = sp.Symbol('b5')
b6 = sp.Symbol('b6')

beta = [1.06, 4.4, 1.2, 11.5, 1.5, 4.0]  # [4.4, 1.2, 11.5, 4.8]
# beta = [1.0629, 540.4003, 1.1596, 11.5964, 11.5964, 4.8254]

# Nonlinearities:
f1 = b6*Wd / ((1 + C)*(1 + b3*Wl)) - Wl

expr_f2 = b1 / (1 + b4*Wl)
f2 = 0.0*(expr_f2 - A)

# expr_f3 = b2/2*Wl + b2/2*S
# expr_f3 = b3*Wl + b4*S
expr_f3 = b2*Wl**2
f3 = expr_f3 - Wd

expr_f4 = Wd / (1 + b5*Wl)
f4 = expr_f4 - C

expr_f5 = Wl
f5 = 0.0*(expr_f5 - S)

f1_eval = f1
for j in range(1,5):
    f1_eval = f1_eval.subs({S: expr_f5, C: expr_f4, Wd: expr_f3, A: expr_f2, b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5]})
u0 = sp.solve(f1_eval, Wl)
print(u0)
u0_real = [np.real(np.complex128(u)) for u in u0 if abs(sp.im(u)) < 10e-15]

if u0_real:
    Wl0 = max(u0_real)
else:
    print("Could not find real steady state!")
    Wl0 = 42 # Some random value

for j in range(1,5):
    expr_f2 = expr_f2.subs({Wd: expr_f3, C: expr_f4, S: expr_f5})
    expr_f3 = expr_f3.subs({A: expr_f2, C: expr_f4, S: expr_f5})
    expr_f4 = expr_f4.subs({A: expr_f2, Wd: expr_f3, S: expr_f5})
    expr_f5 = expr_f5.subs({A: expr_f2, Wd: expr_f3, C: expr_f4})

A0 = expr_f2.subs({Wl: Wl0, b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5]})
Wd0 = expr_f3.subs({Wl: Wl0, b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5]})
C0 = expr_f4.subs({Wl: Wl0, b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5]})
S0 = expr_f5.subs({Wl: Wl0, b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5]})

# test steady state:
print("Testing steady state:")
print(f"Steady state values: Wl={Wl0}, A={A0}, Wd={Wd0}, C={C0}, S={S0}")
print(f1.subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0}))
print(f2.subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0}))
print(f3.subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0}))
print(f4.subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0}))
print(f5.subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0}))
print("")

j11 = sp.diff(f1, Wl).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j12 = sp.diff(f1, A).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j13 = sp.diff(f1, Wd).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j14 = sp.diff(f1, C).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j15 = sp.diff(f1, S).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j21 = sp.diff(f2, Wl).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j22 = sp.diff(f2, A).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j23 = sp.diff(f2, Wd).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j24 = sp.diff(f2, C).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j25 = sp.diff(f2, S).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j31 = sp.diff(f3, Wl).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j32 = sp.diff(f3, A).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j33 = sp.diff(f3, Wd).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j34 = sp.diff(f3, C).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j35 = sp.diff(f3, S).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j41 = sp.diff(f4, Wl).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j42 = sp.diff(f4, A).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j43 = sp.diff(f4, Wd).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j44 = sp.diff(f4, C).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j45 = sp.diff(f4, S).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j51 = sp.diff(f5, Wl).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j52 = sp.diff(f5, A).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j53 = sp.diff(f5, Wd).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j54 = sp.diff(f5, C).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})
j55 = sp.diff(f5, S).subs({b1: beta[0], b2: beta[1], b3: beta[2], b4: beta[3], b5: beta[4], b6: beta[5], Wl: Wl0, Wd: Wd0, C: C0, A: A0, S: S0})

Jacobian = sp.Matrix([[j11, j12, j13, j14, j15],
                      [j21, j22, j23, j24, j25],
                      [j31, j32, j33, j34, j35],
                      [j41, j42, j43, j44, j45],
                      [j51, j52, j53, j54, j55]])
print("Jacobian Matrix at the equilibrium point:")
sp.pprint(Jacobian)
print("")
print("Eigenvalues of the Jacobian Matrix:")
eigenvals = Jacobian.eigenvals()
for val, mult in eigenvals.items():
    print(f"Eigenvalue: {val}, Multiplicity: {mult}")
print("")
print("Jacobian Matrix at the equilibrium point for copying:")
print(Jacobian)
Jac_np = np.array(Jacobian).astype(np.float64)
print('Check Jacobian for unstable subsystems:')
principal_submatrices = get_principle_submatrices(Jac_np)
submatrices_stable = []
for dim, submats in principal_submatrices.items():
    print('')
    print(f"Dimension: {dim}x{dim}")
    print('Stability:')
    submatrices_stable.append(get_stability_submatrices(submats))
    # print('All submatrices:')
    # for submatrix in submats:
    #     print(submatrix)
print('')
if submatrices_stable == [True for dim in principal_submatrices.keys()]:
    print('All principal submatrices are stable!')
else:
    print('There is at least one unstable principal submatrix.')

# Jac_variant_2 = np.array([[-1.54835418090850, -0.968673468738539, 0, -0.519279522931973, 1.00000000000000], [-0.0713698848523722, -1, 0, 0, 0], [1.15960000000000, 0, -1, 0, 11.5964000000000], [-0.896792686154216, 0, 0.0760938750398222, -1, 0], [1, 0, 0, 0, -1]])
