import sympy as sp
from scipy.optimize import fsolve

##################################################################################
# Help functions to work with symbolic functions

def symbolic_to_executable_three_vars(sympy_func):
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
    x, y, z = sp.symbols('x y z')

    # Lambdify the sympy function
    func = sp.lambdify((x, y, z), sympy_func, 'numpy')

    return func


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


def get_constant_steady_state(func_1, func_2, func_3, initial_guess):
    """
    Computes the constant steady state by solving the system of nonlinear equations.
        
    Parameters:
    func_1, func_2, func_3 : string
        The 3 equations as functions of 3 variables (x, y, z).
    initial_guess : tuple
        An initial guess for the solution, (x^0, y^0, z^0).
        
    Returns:
    tuple
        The steady-state solution (x, y, z).
    """
    function_1 = symbolic_to_executable_three_vars(func_1)
    function_2 = symbolic_to_executable_three_vars(func_2)
    function_3 = symbolic_to_executable_three_vars(func_3)
    
    solution = solve_system_of_equations_numerically([function_1, function_2, function_3], initial_guess)
    return solution[0], solution[1], solution[2]


def evaluate_symbolic_at_value_three_vars(sympy_func, x, y, z):
    """
    Evaluates a symbolic sympy function in three variables at given values.

    Parameters:
    sympy_func : sympy expression
        The symbolic function to be evaluated.
    x, y, z : float
        The values at which the function will be evaluated.

    Returns:
    float
        The result of the function evaluation.
    """
    func = symbolic_to_executable_three_vars(sympy_func)
    return func(x, y, z)

