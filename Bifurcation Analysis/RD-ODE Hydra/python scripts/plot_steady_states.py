from scipy.optimize import fsolve
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def solve_system_of_equations_numerically(funcs, initial_guess):
    """
    Solves a system of n equations in n variables numerically.

    Parameters:
    funcs : list of callable
        The n equations as functions of n variables (x_1, ..., x_n).
    initial_guess : tuple
        An initial guess for the solution, e.g., (x^0_1, ..., x^0_n).

    Returns:
    tuple
        A tuple containing the solution (x_1, ..., x_n).
    """
    def equations(variables):
        return [func(*variables) for func in funcs]

    solution = fsolve(equations, initial_guess)
    return solution

def get_constant_steady_state(func_1, func_2, initial_guess):
    """
    Computes the constant steady state by solving the system of nonlinear equations.
        
    Parameters:
    func_1, func_2 : string
        The 2 equations as functions of 2 variables (x, y).
    initial_guess : tuple
        An initial guess for the solution, (x^0, y^0).
        
    Returns:
    tuple
        The steady-state solution (x, y).
    """
    function_1 = symbolic_to_executable_two_vars(func_1)
    function_2 = symbolic_to_executable_two_vars(func_2)
    
    solution = solve_system_of_equations_numerically([function_1, function_2], initial_guess)
    return solution[0], solution[1]

def symbolic_to_executable_two_vars(sympy_func):
    """
    Transforms a symbolic sympy function in three variables to an executable function using lambdify.

    Parameters:
    sympy_func : sympy expression
        The symbolic function to be transformed.

    Returns:
    callable
        A function that can be evaluated numerically.
    """
    # Define the symbols
    x, y = sp.symbols('x y')

    # Lambdify the sympy function
    func = sp.lambdify((x, y), sympy_func, 'numpy')

    return func


if __name__ == '__main__':
    a = 1.5
    b = 2.0
    muu = 0.5
    muv = 1.0
    func_g = f"{b}*x**2 - {muv}*y"
    v0 = (a*muv + np.sqrt(a**2 * muv**2 - 4 * b * muu**2 *muv)) / (2 * b * muu)
    initial_guess_steady_state = [v0, b / muv * v0**2]

    steady_states = []

    for pu in np.arange(0.0, 0.1, 0.0005):
        func_f = f"({a}*x**2 + {pu}) / (1 + y) - {muu}*x"
        
        # Find all steady states for this pu using sympy:
        x, y = sp.symbols('x y')
        sympy_func_f = sp.sympify(func_f)
        sympy_func_g = sp.sympify(func_g)
        sympy_solution = sp.solve([sympy_func_f, sympy_func_g], (x, y))
        print(f"pu: {pu}, sympy steady states: {sympy_solution}")
        steady_states.append(sympy_solution)

    # Plot the steady state x as a function of pu if they are real and positive:
    pu_values = np.arange(0.0, 0.1, 0.0005)
    x_values = []
    y_values = []
    for i, pu in enumerate(pu_values):
        for sol in steady_states[i]:
            if abs(sp.im(sol[0])) < 1e-10 and abs(sp.im(sol[1])) < 1e-10:
                x_values.append(float(sp.re(sol[0])))
                y_values.append(float(sp.re(sol[1])))
                plt.plot(pu, float(sp.re(sol[0])), 'bo')  # Plot x vs pu
    plt.xlabel('pu')
    plt.ylabel('Steady state values for u')
    plt.title('Steady state values for u as a function of pu')
    plt.show()