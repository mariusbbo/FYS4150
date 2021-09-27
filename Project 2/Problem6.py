import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})

"""
Problem 6:
- Read N and number of iterations from file.
- Plot number of iterations against N and estimate by
testing values for variable y.
- Plot number of iterations against N for tridiagonal and dense matrix
"""

# read N and number of iterations from file and create arrays
file = open("Nvals_iterations.txt", "r")
lines = np.array(file.read().split(), dtype=float)
N = lines[::2]
iterations = lines[1::2]

file_rand = open("Nvals_iterations_random.txt", "r")
lines_rand = np.array(file_rand.read().split(), dtype=float)
N_rand = lines_rand[::2]
iterations_rand = lines_rand[1::2]

# test/estimation
y = 1.8*N**2

# plot number of iteerations as function of N and y
plt.plot(N, iterations, label="Tridiagonal matrix")
plt.plot(N, y, label="$1.8N^2$")
plt.title("Number of iterations as function of size of matrix N")
plt.xlabel("N")
plt.ylabel("Iterations")
plt.legend()
plt.show()

# plot number of iterations against N for tridiagonal and dense matrix
plt.plot(N, iterations, label="Tridiagonal matrix")
plt.plot(N, iterations_rand, label="Dense matrix")
plt.title("Number of iterations as function of size of matrix N")
plt.xlabel("N")
plt.ylabel("Iterations")
plt.legend()
plt.show()
