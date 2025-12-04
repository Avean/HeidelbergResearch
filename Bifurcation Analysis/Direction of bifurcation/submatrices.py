import numpy as np
from numpy import linalg as LA
from itertools import combinations

##################################################################################
# Submatrices help functions
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
                print(f'Submatrix\n{submatrix}\nis unstable!')
                all_submatrices_stable = False
                break
    return all_submatrices_stable


def check_for_unstable_submatrix(Jacobian):
    dict_submatrices = get_principle_submatrices(Jacobian)

    submatrices_stable = []
    for dim, submats in dict_submatrices.items():
        # submats list of all submatrices of dimension dim
        submatrices_stable.append(get_stability_submatrices(submats))
    if submatrices_stable == [True for _ in dict_submatrices.keys()]:
        print('All principal submatrices are stable!')
        return False
    else:
        print('There is at least one unstable principal submatrix.')
        return True

