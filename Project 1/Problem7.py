import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})

file_u = open("x_u_vals.txt", "r")
lines_u = np.array(file_u.read().split(), dtype=float)
x = lines_u[::2]
u = lines_u[1::2]

file_v1 = open("x_v_vals_n10.txt", "r")
lines_v1 = np.array(file_v1.read().split(), dtype=float)
x_v1 = lines_v1[::2]
v1 = lines_v1[1::2]

file_v2 = open("x_v_vals_n100.txt", "r")
lines_v2 = np.array(file_v2.read().split(), dtype=float)
x_v2 = lines_v2[::2]
v2 = lines_v2[1::2]

file_v3 = open("x_v_vals_n1000.txt", "r")
lines_v3 = np.array(file_v3.read().split(), dtype=float)
x_v3 = lines_v3[::2]
v3 = lines_v3[1::2]

file_v4 = open("x_v_vals_n10000.txt", "r")
lines_v4 = np.array(file_v4.read().split(), dtype=float)
x_v4 = lines_v4[::2]
v4 = lines_v4[1::2]

# plot comparing v with u for different n
plt.plot(x, u, label="u(x)")
plt.plot(x_v1, v1, label="v(x), n=10")
plt.plot(x_v2, v2, label="v(x), $n=10^2$")
plt.plot(x_v3, v3, label="v(x), $n=10^3$")
plt.plot(x_v4, v4, label="v(x), $n=10^4$")
plt.title("u(x) and v(x) plotted for different n")
plt.xlabel("x")
plt.ylabel("u(x) (v(x))")
plt.legend()
#plt.savefig("comp_u_v_different_n.pdf")
plt.show()
