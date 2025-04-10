import numpy as np
import matplotlib.pyplot as plt
import time

N = 2000
A = np.random.rand(N,N)
B = np.random.rand(N,N)


def fun(A,B):
    return np.matmul(A,B)


t1 = time.perf_counter()
C = fun(A,B)
t2 = time.perf_counter()

print("Time taken for matrix multiplication: ", t2-t1)


