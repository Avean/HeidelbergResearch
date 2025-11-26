import numpy as np
from itertools import combinations
from numpy import linalg as LA

A = np.array([
    [-0.001, -0.0007, -0.0007, -0.0028],
    [0.004, -0.001, -0.0007, -0.0028],
    [0.004, 0.004, -0.001, -0.0028],
    [-1.0, -1.0, -1.0, -4.0]
])

# A = np.array([
#     [-0.7, -0.0007, -0.0007, 0.1],
#     [0.0007, -0.7, 0.1, -0.01],
#     [-2., -4., -8., -1.],
#     [1., 1., 1., -4.]
# ])

Jac = np.array([
    [-1.08652, -0.0528344, 0.0,      -0.0286439, 1.0],
    [-3.25103, -1.0,       0.0,      0.0,        0.0],
    [44.1375,  0.0,        -1.0,     0.0,        44.1375],
    [-11.0262, 0.0,        0.513573, -1.0,       0.0],
    [1.0,      0.0,        0.0,      0.0,        -1.0]
])

# Jac_variant = np.array([[-1.0030148996730284, -0.0012835606066773457, 0.0, -0.0011013202189362518, 1.000000000000001],
#                [-11.61282169573653, -1.0, 0.0, 0.0, 0.0],
#                [270.20015, 0.0, -1.0, 0.0, 270.20015],
#                [-15.397028807704197, 0.0, 0.9706464752082521, -1.0, 0.0],
#                [1.0, 0.0, 0.0, 0.0, -1.0]])
# Jac_variant_2 = np.array([[-1.54835418090850, -0.968673468738539, 0, -0.519279522931973, 1.00000000000000], [-0.0713698848523722, -1, 0, 0, 0], [1.15960000000000, 0, -1, 0, 11.5964000000000], [-0.896792686154216, 0, 0.0760938750398222, -1, 0], [1, 0, 0, 0, -1]])
# Jac_variant_2 = np.array([[-1.54835418090850, -0.968673468738539, 0, -0.519279522931973, 1.00000000000000], [-0.0713698848523722, -1, 0, 0, 0], [1.15960000000000, 0, -1, 0, 11.5964000000000], [-0.896792686154216, 0, 0.0760938750398222, -1, 0], [1, 0, 0, 0, -1]])

# Jac_variant_2 = np.array([[-1.00301489967303, -0.00128356060667735, 0, -0.00110132021893625, 1.00000000000000], [-11.6128216957365, -1, 0, 0, 0], [540.400300000000, 0, -1, 0, 0], [-15.3970288077042, 0, 0.970646475208252, -1, 0], [1, 0, 0, 0, -1]])
Jac_variant_2 = np.array([[-1.00301489967303, -0.00128356060667735, 0.00185048009780898, -0.00110132021893625, 0], [-11.6128216957365, -1, 0, 0, 0], [540.400300000000, 0, -1, 0, 0], [-15.3970288077042, 0, 0.970646475208252, -1, 0], [1, 0, 0, 0, -1]])
Jac_variant_2 = np.array([[-1.00301489967303, -0.00128356060667735, 0.00185048009780898, -0.00110132021893625, 0], [-11.6128216957365, -1, 0, 0, 0], [540.400300000000, 0, -1, 0, 0], [-15.3970288077042, 0, 0.970646475208252, -1, 0], [0, 0, 0, 0, 0]])
# Kickout A, replace third equation with linear term b2*Wloc - Wdif and replace S in first equation with Wdif/b2 and kickout S:
Jac_variant_2 = np.array([[-1.00875780675000, 0.00185048009780898, -0.00159292112319610], [540.400300000000, -1, 0], [-40.3092835636705, 0.918817947333312, -1]])
# Kickout A and replace third equation with linear term b2*Wloc - Wdif
Jac_variant_2 = np.array([[-1.00875780675000, 0, -0.00159292112319610, 1.00000000000000], [540.400300000000, -1, 0, 0], [-40.3092835636705, 0.918817947333312, -1, 0], [1, 0, 0, -1]])
def get_submatrices(matrix):
    n = len(matrix)
    submatrices = {}
    for dim in range(1, n + 1):
        submatrices[dim] = []
        for i in range(n - dim + 1):
            for j in range(n - dim + 1):
                submatrix = matrix[i:i + dim, j:j + dim]
                submatrices[dim].append(submatrix)
    return submatrices


def get_principle_submatrices(matrix):
    n = len(matrix)
    submatrices = dict()
    for dim in range(1, n + 1):
        submatrices[dim] = []
        for indices in combinations(range(n), dim):
            submatrix = matrix[np.ix_(indices, indices)]
            submatrices[dim].append(submatrix)

    return submatrices


def get_stability_submatrices(submatrices):
    all_submatrices_stable = True
    for submatrix in submatrices:
        eigvs, _ = LA.eig(submatrix)
        for eig in eigvs:
            if np.real(eig) > 0:
                print('Submatrix ' + str(submatrix) + ' is unstable!')
                all_submatrices_stable = False
                break
    if all_submatrices_stable:
        print('Submatrices from this list not unstable.')
    return all_submatrices_stable


if __name__ == '__main__':
    principal_submatrices = get_principle_submatrices(Jac_variant_2)
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

    # print('Unstable are J125, J145, J1235, J1245')
    # J125 = Jac_variant_2[np.ix_([0,1,4], [0,1,4])]
    # print(J125)
    # print(LA.eig(J125)[0])
    # print('')
    # J145 = Jac_variant_2[np.ix_([0,3,4], [0,3,4])]
    # print(J145)
    # print(LA.eig(J145)[0])
    # print('')
    # J1235 = Jac_variant_2[np.ix_([0,1,2,4], [0,1,2,4])]
    # print(J1235)
    # print(LA.eig(J1235)[0])
    # print('')
    # J1245 = Jac_variant_2[np.ix_([0,1,3,4], [0,1,3,4])]
    # print(J1245)
    # print(LA.eig(J1245)[0])
    # print('')


    # print('')
    # print('')
    # print('Other submatrices (not principal):')
    # submatrices = get_submatrices(A)
    # for dim, submats in submatrices.items():
    #     print('')
    #     print(f"Dimension: {dim}x{dim}")
    #     print('Stability:')
    #     get_stability_submatrices(submats)
    #     print('All submatrices:')
    #     for submatrix in submats:
    #         print(submatrix)
