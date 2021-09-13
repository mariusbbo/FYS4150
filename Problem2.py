import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})

file_u = open("x_u_vals.txt", "r")
lines_u = np.array(file_u.read().split(), dtype=float)
x = lines_u[::2]
u = lines_u[1::2]

# plot exact solution u(x) as function of x
plt.plot(x, u, label="u(x)")
plt.title("Plot of exact solution, u(x)")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.legend()
#plt.savefig("u(x).pdf")
plt.show()
