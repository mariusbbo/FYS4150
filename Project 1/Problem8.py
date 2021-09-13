import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})

file_u1 = open("x_u_vals_n10.txt", "r")
lines_u1 = np.array(file_u1.read().split(), dtype=float)
x_u1 = lines_u1[::2]
u1 = lines_u1[1::2]

file_u2 = open("x_u_vals_n100.txt", "r")
lines_u2 = np.array(file_u2.read().split(), dtype=float)
x_u2 = lines_u2[::2]
u2 = lines_u2[1::2]

file_u3 = open("x_u_vals_n1000.txt", "r")
lines_u3 = np.array(file_u3.read().split(), dtype=float)
x_u3 = lines_u3[::2]
u3 = lines_u3[1::2]

file_u4 = open("x_u_vals_n10000.txt", "r")
lines_u4 = np.array(file_u4.read().split(), dtype=float)
x_u4 = lines_u4[::2]
u4 = lines_u4[1::2]

file_u5 = open("x_u_vals_n10^5.txt", "r")
lines_u5 = np.array(file_u5.read().split(), dtype=float)
x_u5 = lines_u5[::2]
u5 = lines_u5[1::2]

file_u6 = open("x_u_vals_n10^6.txt", "r")
lines_u6 = np.array(file_u6.read().split(), dtype=float)
x_u6 = lines_u6[::2]
u6 = lines_u6[1::2]

file_u7 = open("x_u_vals_n10^7.txt", "r")
lines_u7 = np.array(file_u7.read().split(), dtype=float)
x_u7 = lines_u7[::2]
u7 = lines_u7[1::2]

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

file_v5 = open("x_v_vals_n10^5.txt", "r")
lines_v5 = np.array(file_v5.read().split(), dtype=float)
x_v5 = lines_v5[::2]
v5 = lines_v5[1::2]

file_v6 = open("x_v_vals_n10^6.txt", "r")
lines_v6 = np.array(file_v6.read().split(), dtype=float)
x_v6 = lines_v6[::2]
v6 = lines_v6[1::2]

file_v7 = open("x_v_vals_n10^7.txt", "r")
lines_v7 = np.array(file_v7.read().split(), dtype=float)
x_v7 = lines_v7[::2]
v7 = lines_v7[1::2]

# absolute error
delta_n10 = abs(u1-v1)
delta_n10_2 = abs(u2-v2)
delta_n10_3 = abs(u3-v3)
delta_n10_4 = abs(u4-v4)
delta_n10_5 = abs(u5-v5)
delta_n10_6 = abs(u6-v6)
delta_n10_7 = abs(u7-v7)

#relative error
# skip first value to avoid division by zero
# skip last value to avoid computational errors not giving zero for
#u(1) meaning the relative error becomes very large when it should be zero
epsilon_n10 = abs(delta_n10[1:-1]/u1[1:-1])
epsilon_n10_2 = abs(delta_n10_2[1:-1]/u2[1:-1])
epsilon_n10_3 = abs(delta_n10_3[1:-1]/u3[1:-1])
epsilon_n10_4 = abs(delta_n10_4[1:-1]/u4[1:-1])
epsilon_n10_5 = abs(delta_n10_5[1:-1]/u5[1:-1])
epsilon_n10_6 = abs(delta_n10_6[1:-1]/u6[1:-1])
epsilon_n10_7 = abs(delta_n10_7[1:-1]/u7[1:-1])


# plot absolute error
# skip first and last value of delta to avoid log10 of zero
plt.plot(x_u1[1:-1], np.log10(delta_n10[1:-1]), label="n=10")
plt.plot(x_u2[1:-1], np.log10(delta_n10_2[1:-1]), label="$n=10^2$")
plt.plot(x_u3[1:-1], np.log10(delta_n10_3[1:-1]), label="$n=10^3$")
plt.plot(x_u4[1:-1], np.log10(delta_n10_4[1:-1]), label="$n=10^4$")
plt.title("Absolute error |$u_i - v_i$| for different n")
plt.xlabel("x")
plt.ylabel("Absolute error")
plt.legend()
plt.savefig("Absolute_error.pdf")
plt.show()

# plot relative error
plt.plot(x_u1[1:-1], np.log10(epsilon_n10), label="n=10")
plt.plot(x_u2[1:-1], np.log10(epsilon_n10_2), label="$n=10^2$")
plt.plot(x_u3[1:-1], np.log10(epsilon_n10_3), label="$n=10^3$")
plt.plot(x_u4[1:-1], np.log10(epsilon_n10_4), label="$n=10^4$")
plt.title(r"Relative error |$\frac{u_i - v_i}{u_i}$| for different n")
plt.xlabel("x")
plt.ylabel("Relative error")
plt.legend()
plt.savefig("Relative_error.pdf")
plt.show()

# list with max relative errors
max_eps_list = [max(epsilon_n10), max(epsilon_n10_2), max(epsilon_n10_3), max(epsilon_n10_4),
max(epsilon_n10_5), max(epsilon_n10_6), max(epsilon_n10_7)]

n_list = [10, 100, 1000, 10**4, 10**5, 10**6, 10**7] # list with n values

# print max relative errors
for i in range(len(n_list)):
    print("Max relative error for {:.0e} is {:1.3e}".format(n_list[i], max_eps_list[i]))

# plot max of relative error as function of n
plt.plot(n_list, max_eps_list)
plt.title("Max relative error for different n")
plt.xlabel("n")
plt.ylabel("Max relative error")
plt.xscale("log")
plt.yscale("log")
plt.savefig("max_relative_error.pdf")
plt.show()
