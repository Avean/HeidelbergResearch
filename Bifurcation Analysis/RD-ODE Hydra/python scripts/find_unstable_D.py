import numpy as np
import matplotlib.pyplot as plt

# Define the matrix A
# A = np.array([
#     [-7.e-01, -7.e-04, -7.e-04,  1.e-01],
#     [ 7.e-04, -7.e-01,  1.e-01, -1.e-02],
#     [-2.e+00, -4.e+00, -8.e+00, -1.e+00],
#     [ 1.e+00,  1.e+00,  1.e+00, -4.e+00]
# ])
A = np.array([
    [-1.e-03, -7.e-04, -7.e-04, -0.0028],
    [ 4.e-03, -1.e-03, -7.e-04, -0.0028],
    [ 4.e-03,  4.e-03, -1.e-03, -0.0028],
    [-1.e+00, -1.e+00, -1.e+00, -4.e+00]
])
Jac = np.array([
    [-1.08652, -0.0528344, 0.0,      -0.0286439, 1.0],
    [-3.25103, -1.0,       0.0,      0.0,        0.0],
    [44.1375,  0.0,        -1.0,     0.0,        44.1375],
    [-11.0262, 0.0,        0.513573, -1.0,       0.0],
    [1.0,      0.0,        0.0,      0.0,        -1.0]
])

Jac_woS = np.array([
        [-0.0865168027146120, -0.0528343931336280, 0, -0.0286439413830729],
        [-3.25102684618881, -1, 0, 0],
        [88.2749069711968, 0, -1, 0],
        [-11.0262316692890, 0, 0.513573351438324, -1]
        ])

Jac_woAS = np.array([[-0.121780213632072, 0, -0.0282181650838908],
                     [129.244146954899, -1, 0],
                     [-15.7313780996255, 0.418985210160453, -1]])

Jac_variant = np.array([[-1.0030148996730284, -0.0012835606066773457, 0.0, -0.0011013202189362518, 1.000000000000001],
               [-11.61282169573653, -1.0, 0.0, 0.0, 0.0],
               [270.20015, 0.0, -1.0, 0.0, 270.20015],
               [-15.397028807704197, 0.0, 0.9706464752082521, -1.0, 0.0],
               [1.0, 0.0, 0.0, 0.0, -1.0]])

Jac_variant = np.array([[-1.00301489967303, -0.00128356060667735, 0, -0.00110132021893625, 1.00000000000000], [-11.6128216957365, -1, 0, 0, 0], [540.400300000000, 0, -1, 0, 0], [-15.3970288077042, 0, 0.970646475208252, -1, 0], [1, 0, 0, 0, -1]])

# Jac_variant_2 = np.array([[-1.00875780675000, 0, -0.00159292112319610, 1.00000000000000], [540.400300000000, -1, 0, 0], [-40.3092835636705, 0.918817947333312, -1, 0], [1, 0, 0, -1]])

# Function to compute eigenvalues for A - D with given d
def compute_eigenvalues(d, k=1):
    # nu = [0.0, 3.8154e-05, 0.4433, 6.0713e-08, 0.0004]
    D = np.diag([0.0, 3.8154e-05, d, 6.0713e-08, 0.0004])
    # D = np.diag([0.0, d, 1e-5, 10.0, 1e-6])
    # D = np.diag([0.0, 3.8154e-03, d, 6.0713e-03])
    # D = np.diag([0.0, d, 6.0713e-08])
    return np.linalg.eigvals(Jac_variant - k**2*np.pi**2*D)


# Function to perform binary search for the critical d value
def find_critical_d(d_min, d_max, tolerance=1e-6, k=1, fast=True):
    if fast:
        while (d_max - d_min) > tolerance:
            d = d_min + d_max*tolerance
            eigenvalues = compute_eigenvalues(d, k)
            if np.any(np.real(eigenvalues) >= 0):
                return d
            else:
                d_min = d
        return d_min
    else:
        while (d_max - d_min) > tolerance:
            d = d_max - tolerance
            eigenvalues = compute_eigenvalues(d, k)
            if np.any(np.real(eigenvalues) >= 0):
                return d
            else:
                d_max = d
        return d_min
    # while (d_max - d_min) > tolerance:
    #     d = (d_max + d_min) / 2
    #     eigenvalues = compute_eigenvalues(d)
    #     if np.any(np.real(eigenvalues) >= 0):
    #         d_max = d
    #     else:
    #         d_min = d
    # return (d_max + d_min) / 2


def find_max_d(d_min, d_max, tolerance=1e-8):
    d_min = find_critical_d(d_min, d_max, tolerance)
    if (d_max - d_min) <= tolerance:
        return d_min
    else:
        while (d_max - d_min) > tolerance:
            d = d_min + d_max * tolerance
            eigenvalues = compute_eigenvalues(d)
            if np.all(np.real(eigenvalues) < 0):
                return d
            else:
                d_min = d
        return d_min


def calculate_eigvalues_for_int(d_start, d_stop, d_num):
    d_values = np.linspace(d_start, d_stop, d_num)
    eigenvalues = np.array([compute_eigenvalues(d) for d in d_values])
    return eigenvalues


if __name__ == '__main__':
    # Narrow the search interval for more precision
    # critical_d = find_critical_d(0, 10, k=8)
    # print(f"The critical value of d is approximately: {critical_d}")
    # max_d = find_max_d(0, 150)
    # print(f"The upper critical value of d is approximately: {max_d}")
    # print(compute_eigenvalues(0))

    # For A:
    # print(compute_eigenvalues(0.12))
    # print(compute_eigenvalues(0.125))
    # print(compute_eigenvalues(1.82))
    # print(compute_eigenvalues(1.83))

    # For Jac:
    # print(compute_eigenvalues(0.0078, k=8))
    # print(compute_eigenvalues(0.008, k=8))
    # print(compute_eigenvalues(0.0088, k=8))
    # print(compute_eigenvalues(0.0094, k=8))

    # print(compute_eigenvalues(0.18, k=5))
    # print(compute_eigenvalues(0.19, k=5))
    # print(compute_eigenvalues(0.20, k=5))
    # print('')
    # print(compute_eigenvalues(0, k=1))
    # print(compute_eigenvalues(0, k=2))
    # print(compute_eigenvalues(0.2, k=1))
    # print(compute_eigenvalues(0.0004, k=24))
    # print(compute_eigenvalues(0.0005, k=24))

    for k in range(1, 13):
        # print(compute_eigenvalues(0.4433, k))
        print(f"k={k}: {find_critical_d(0, 2.3, tolerance=1e-4, k=k, fast=True)}")

    # # print(compute_eigenvalues(4.6e-07, 1440))
    # # print(compute_eigenvalues(4.7e-07, 1440))
    # # print(compute_eigenvalues(4.8e-07, 1450))
    # # print(f"k={1440}: {find_critical_d(0, 0.5, tolerance=1e-9, k=1440, fast=True)}")
    # # print(f"k={1450}: {find_critical_d(0, 0.5, tolerance=1e-9, k=1450, fast=True)}")
    # # print(f"k={1460}: {find_critical_d(0, 0.5, tolerance=1e-9, k=1455, fast=True)}")
    # # D = np.diag([0.0, 3.8154e-05, 0.4433, 6.0713e-08, 0.0])
    # D = np.diag([0.0, 3.8154e-03, 0.4433, 6.0713e-03])
    # # D = np.diag([0.0, 0.4433, 6.0713e-08])
    # # print(np.linalg.det(Jac_woS - 2104**2*np.pi**2*D))
    # # print(np.linalg.det(Jac_woS - 2105**2*np.pi**2*D))
    # list_k = np.linspace(0,13,100)
    # list_det = [-np.linalg.det(Jac_woS - k**2*np.pi**2*D) for k in list_k]
    # A = [x for x in list_det if x > 0]
    # list_k_pos = [k for (ind, k) in enumerate(list_k) if list_det[ind] > 0]
    # # plt.plot(list_k, list_det)
    # plt.plot(list_k_pos, A)
    # plt.show()
    # # print(list_det)