import sympy as sp
from sympy import Matrix

def get_bif_points_fast(J, D, max_wn = 200):
    # d = sp.Symbol("d")
    wn = sp.Symbol("wn")
    Mw = J - wn**2 * sp.pi**2 * D
    detw = Mw.det()
    dw = sp.solve(detw, d)[0]
    # print(dw)
    bif_points = []
    wave_numbers = []
    for k in range(1, max_wn + 1):
        dk = dw.subs(wn, k)
        if dk > 0:
            bif_points.append(dk)
            wave_numbers.append(k)
        else:
            break
    return wave_numbers, bif_points


if __name__ == '__main__':

    # Jac = Matrix([
    #     [-1.0865168027146117, -0.05283439313362805, 0.0, -0.028643941383072864, 0.9999999999999997],
    #     [-3.2510268461888066, -1.0, 0.0, 0.0, 0.0],
    #     [44.1374534855984, 0.0, -1.0, 0.0, 44.1374534855984],
    #     [-11.026231669288988, 0.0, 0.5135733514383242, -1.0, 0.0],
    #     [1.0, 0.0, 0.0, 0.0, -1.0]
    #     ])
    d = sp.Symbol("d")
    # D = Matrix.diag(0.0, 3.8154e-05, d, 6.0713e-08, 0.0004)
    # D = Matrix.diag(0.0, 0.006980161761761762, d, 0.0010765240580580579, 0.004522162162162163)


    # Jac = Matrix([[-1.98080925169869, 0.00533624701923220, -2.46725578935589], [374.795243322107, -1, 0], [-0.381051127243924, 0.00203753448872319, -1]])
    # D = Matrix.diag(0.0, d, 6.5)

    Jac = Matrix([[-1.84228031377574, 0.0510690552489811, -0.360247875064841], [39.1626590750355, -1, 0], [-2.21878322298442, 0.130285465527934, -1]])
    Jac = Matrix([[-1.81327685161692, 0.183674868415938, -0.892587063306863], [10.8888059496029, -1, 0], [-0.713731963241883, 0.155173412325415, -1]])
    D = Matrix.diag(0.0, d, 1.0)

    wave_numbers, bif_points = get_bif_points_fast(Jac, D, max_wn=5)
    print(f"{len(bif_points)} bifurcation points {bif_points} cor. to wave numbers: {wave_numbers}")