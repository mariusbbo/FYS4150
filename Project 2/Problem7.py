import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})

"""
Problem 7
- Read files and plot eigenvectors
for n = 10 and n = 100
"""
# read numerical eigenvectors from file and create arrays
file = open("numerical_eigenvectors_n10.txt", "r")
lines = np.array(file.read().split(), dtype=float)
v1_ = lines[::3]
v2_ = lines[1::3]
v3_ = lines[2::3]*-1
# scale eigenvector v3 to get same sign on elements as in analytical eigenvector

# add boundary conditions
v1_n10v = np.concatenate(([0], v1_, [0]))
v2_n10v = np.concatenate(([0], v2_, [0]))
v3_n10v = np.concatenate(([0], v3_, [0]))

# read analytical eigenvectors from file and create arrays
file = open("analytical_eigenvectors_n10.txt", "r")
lines = np.array(file.read().split(), dtype=float)
u1_ = lines[::3]
u2_ = lines[1::3]
u3_ = lines[2::3]

# add boundary conditions
u1_n10u = np.concatenate(([0], u1_, [0]))
u2_n10u = np.concatenate(([0], u2_, [0]))
u3_n10u = np.concatenate(([0], u3_, [0]))

# read numerical eigenvectors from file and create arrays
file = open("numerical_eigenvectors_n100.txt", "r")
lines = np.array(file.read().split(), dtype=float)
v1_ = lines[::3]
v2_ = lines[1::3]
v3_ = lines[2::3]

# add boundary conditions
v1_n100v = np.concatenate(([0], v1_, [0]))
v2_n100v = np.concatenate(([0], v2_, [0]))
v3_n100v = np.concatenate(([0], v3_, [0]))

# read analytical eigenvectors from file and create arrays
file = open("analytical_eigenvectors_n100.txt", "r")
lines = np.array(file.read().split(), dtype=float)
u1_ = lines[::3]
u2_ = lines[1::3]
u3_ = lines[2::3]

# add boundary conditions
u1_n100u = np.concatenate(([0], u1_, [0]))
u2_n100u = np.concatenate(([0], u2_, [0]))
u3_n100u = np.concatenate(([0], u3_, [0]))

# array with values for dimensionless length scale x^hat
x_n10 = np.linspace(0, 1, 11)
x_n100 = np.linspace(0, 1, 101)

# plot to compare numerical and analytical eigenvectors

fig,ax = plt.subplots(1,2)

ax[0].plot(x_n10, u1_n10u, label="u1")
ax[0].plot(x_n10, v1_n10v, label="v1")
ax[0].set_title(r"$\vec{v}_1$, $\vec{v}_1$, n = 10")
ax[0].set_xlabel("$\hat{x}$")
ax[0].set_ylabel("u / v")
ax[0].legend()

ax[1].plot(x_n100, u1_n100u, label="u1")
ax[1].plot(x_n100, v1_n100v, label="v1")
ax[1].set_title(r"$\vec{v}_1$, $\vec{u}_1$, n = 100")
ax[1].set_xlabel("$\hat{x}$")
ax[1].legend()

plt.show()

fig, ax = plt.subplots(1,2)

ax[0].plot(x_n10, u2_n10u, label="u2")
ax[0].plot(x_n10, v2_n10v, label="v2")
ax[0].set_title(r"$\vec{v}_2$, $\vec{u}_2$ , n = 10")
ax[0].set_xlabel("$\hat{x}$")
ax[0].set_ylabel("u / v")
ax[0].legend()

ax[1].plot(x_n100, u2_n100u, label="u2")
ax[1].plot(x_n100, v2_n100v, label="v2")
ax[1].set_title(r"$\vec{v}_2$, $\vec{u}_2$, n = 100")
ax[1].set_xlabel("$\hat{x}$")
ax[1].legend()

plt.show()

fig, ax = plt.subplots(1,2)

ax[0].plot(x_n10, u3_n10u, label="u3")
ax[0].plot(x_n10, v3_n10v, label="v3")
ax[0].set_title(r"$\vec{v}_3$, $\vec{u}_3$, n = 10")
ax[0].set_xlabel("$\hat{x}$")
ax[0].set_ylabel("u / v")
ax[0].legend()

ax[1].plot(x_n100, u3_n100u, label="u3")
ax[1].plot(x_n100, v3_n100v, label="v3")
ax[1].set_title(r"$\vec{v}_3$, $\vec{u}_3$, n = 100")
ax[1].set_xlabel("$\hat{x}$")
ax[1].legend(loc="upper right")

plt.show()
