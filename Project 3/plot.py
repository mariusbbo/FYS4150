import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa
from mpl_toolkits.mplot3d import Axes3D
plt.rcParams.update({'font.size': 12})

# function loading data and creating array if the object in
# the file is a cube (3-dimensional array)
def load_file_cube(filename):
    object = pa.cube()
    object.load("{}".format(filename))
    return np.array(object)

# function loading data and creating array if the object in
# the file is a matrix (2-dimensional array)
def load_file_mat(filename):
    object = pa.mat()
    object.load("{}".format(filename))
    return np.array(object)

# plotting function
def plot(x, y, legend, title, xlabel, ylabel):
    plt.plot(x, y, label="{}".format(legend))
    plt.xlabel("{}".format(xlabel))
    plt.ylabel("{}".format(ylabel))
    plt.title("{}".format(title))

# load data and create arrays with positions and velocities of two particles without interactions
pos = load_file_cube("pos.dat")
vel = load_file_cube("vel.dat")
pos_with_ints = load_file_cube("pos_two_particles_with_ints.dat")
vel_with_ints = load_file_cube("vel_two_particles_with_ints.dat")

t = np.linspace(0, 100, len(pos[:,0,0])) # time list

# plot of motion in z-axis as function of time
plot(t, pos[:,2,0], None, r"z-position as function of time", r"Time [$\mu$s]", r"z [$\mu$m]")
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("z_pos_t.pdf")
# plt.show()

# compare the angular frequency from the plot with
# omega_z found from the expression in the project description
max = np.where(pos[10:200000,2,0] > 1-1e-9) # find indices of two peaks in z-time plot
delta_t = t[max[0][1]] - t[max[0][0]] # difference between time of two peaks
freq = 2*np.pi/delta_t # compute angular frequency from z-time plot
print("Angular frequency found from plot: {:.4f}".format(freq))

# particle paramters and Penning trap properties
q = 1
m = 40
B0 = 9.65*10
V0 = 9.65*10**8
d = 10**4
omega_z = np.sqrt(2*q*V0/(m*d**2))
print("omega_z: {:.4f}".format(omega_z))

# plot of motion of two particles in xy-plane without particle interactions
plot(pos[:,0,0], pos[:,1,0], "Particle 1", "Without particle interactions", r"x [$\mu$m]", r"y [$\mu$m]")
plot(pos[:,0,1], pos[:,1,1], "Particle 2", "Without particle interactions", r"x [$\mu$m]", r"y [$\mu$m]")
plt.legend()
plt.axis("equal")
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("two_particles_xy_noints.pdf")
# plt.show()

# plot of motion of two particles in xy-plane including particle interactions
plot(pos_with_ints[:,0,0], pos_with_ints[:,1,0], "Particle 1", "With particle interactions", r"x [$\mu$m]", r"y [$\mu$m]")
plot(pos_with_ints[:,0,1], pos_with_ints[:,1,1], "Particle 2", "With particle interactions", r"x [$\mu$m]", r"y [$\mu$m]")
plt.legend()
plt.axis("equal")
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("two_particles_xy_with_ints.pdf")
# plt.show()

# phase space plots of two particles - x - without interactions
plot(pos[:,0,0], vel[:,0,0], "Particle 1", "Without particle interactions", r"x [$\mu$m]", r"$v_x$ [$\mu$m / $\mu$s]")
plot(pos[:,0,1], vel[:,0,1], "Particle 2", "Without particle interactions", r"x [$\mu$m]", r"$v_x$ [$\mu$m / $\mu$s]")
plt.legend(loc="upper right")
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("phase_space_x_noints.pdf")
# plt.show()

# phase space plots of two particles - x - including interactions
plot(pos_with_ints[:,0,0], vel_with_ints[:,0,0], "Particle 1", "With particle interactions", r"x [$\mu$m]", r"$v_x$ [$\mu$m / $\mu$s]")
plot(pos_with_ints[:,0,1], vel_with_ints[:,0,1], "Particle 2", "With particle interactions", r"x [$\mu$m]", r"$v_x$ [$\mu$m / $\mu$s]")
plt.legend(loc="upper right")
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("phase_space_x_with_ints.pdf")
# plt.show()

# phase space plots of two particles - y - without interactions
plot(pos[:,1,0], vel[:,1,0], "Particle 1", "Without particle interactions", r"y [$\mu$m]", r"$v_y$ [$\mu$m / $\mu$s]")
plot(pos[:,1,1], vel[:,1,1], "Particle 2", "Without particle interactions", r"y [$\mu$m]", r"$v_y$ [$\mu$m / $\mu$s]")
plt.legend(loc="upper right")
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("phase_space_y_noints.pdf")
# plt.show()

# phase space plots of two particles - y - including interactions
plot(pos_with_ints[:,1,0], vel_with_ints[:,1,0], "Particle 1", "With particle interactions", r"y [$\mu$m]", r"$v_y$ [$\mu$m / $\mu$s]")
plot(pos_with_ints[:,1,1], vel_with_ints[:,1,1], "Particle 2", "With particle interactions", r"y [$\mu$m]", r"$v_y$ [$\mu$m / $\mu$s]")
plt.legend(loc="upper right")
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("phase_space_y_with_ints.pdf")
# plt.show()

# phase space plots of two particles - z - without interactions
plot(pos[:,2,0], vel[:,2,0], "Particle 1", "Without particle interactions", r"z [$\mu$m]", r"$v_z$ [$\mu$m / $\mu$s]")
plot(pos[:,2,1], vel[:,2,1], "Particle 2", "Without particle interactions", r"z [$\mu$m]", r"$v_z$ [$\mu$m / $\mu$s]")
plt.legend(loc="upper right")
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("phase_space_z_noints.pdf")
# plt.show()

# phase space plots of two particles - z - including interactions
plot(pos_with_ints[:,2,0], vel_with_ints[:,2,0], "Particle 1", "With particle interactions", r"z [$\mu$m]", r"$v_z$ [$\mu$m / $\mu$s]")
plot(pos_with_ints[:,2,1], vel_with_ints[:,2,1], "Particle 2", "With particle interactions", r"z [$\mu$m]", r"$v_z$ [$\mu$m / $\mu$s]")
plt.legend(loc="upper right")
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("phase_space_z_with_ints.pdf")
# plt.show()

# 3D plot of trajectories of two particles without particle interactions
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot(pos[:,0,0], pos[:,1,0], pos[:,2,0], label="Particle 1")
ax.plot(pos[:,0,1], pos[:,1,1], pos[:,2,1], label="Particle 2")
ax.set_title("Without particle interactions")
ax.set_xlabel(r"x [$\mu$m]")
ax.set_ylabel(r"y [$\mu$m]")
ax.set_zlabel(r"z [$\mu$m]")
ax.xaxis.labelpad=10
ax.yaxis.labelpad=10
ax.zaxis.labelpad=10
plt.legend()
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("two_particles_3D_noints.pdf")
# plt.show()

# 3D plot of trajectories of two particles including particle interactions
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot(pos_with_ints[:,0,0], pos_with_ints[:,1,0], pos_with_ints[:,2,0], label="Particle 1")
ax.plot(pos_with_ints[:,0,1], pos_with_ints[:,1,1], pos_with_ints[:,2,1], label="Particle 2")
ax.set_title("With particle interactions")
ax.set_xlabel(r"x [$\mu$m]")
ax.set_ylabel(r"y [$\mu$m]")
ax.set_zlabel(r"z [$\mu$m]")
ax.xaxis.labelpad=10
ax.yaxis.labelpad=10
ax.zaxis.labelpad=10
plt.legend()
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("two_particles_3D_with_ints.pdf")
# plt.show()

# load position data from RK4
pos_dt_1 = load_file_cube("pos_dt_1.dat")
pos_dt_01 = load_file_cube("pos_dt_0.1.dat")
pos_dt_001 = load_file_cube("pos_dt_0.01.dat")
pos_dt_0001 = load_file_cube("pos_dt_0.001.dat")
pos_dt_00001 = load_file_cube("pos_dt_0.0001.dat")

# # load position data from forward Euler
pos_dt_1_Euler = load_file_cube("pos_dt_1_Euler.dat")
pos_dt_01_Euler = load_file_cube("pos_dt_0.1_Euler.dat")
pos_dt_001_Euler = load_file_cube("pos_dt_0.01_Euler.dat")
pos_dt_0001_Euler = load_file_cube("pos_dt_0.001_Euler.dat")
pos_dt_00001_Euler = load_file_cube("pos_dt_0.0001_Euler.dat")

# function computing relative error and maximum of absolute error (max(delta))
def relative_error(r_num, r_exact):
    t = np.linspace(0, 100, len(r_exact[0])) # time list
    delta = np.zeros(len(r_exact[0])) # absolute error (delta) list
    rel_error = np.zeros(len(r_exact[0])) # relative error list

    # compute absolute and relative error for all positions
    # transpose is used for absolute error (delta) because of shape of analytical solution
    for i in range(len(r_exact[0])):
        delta[i] = np.linalg.norm(np.transpose(r_num[i])[0] - r_exact[:,i]) # compute absolute error
        rel_error[i] = np.abs(delta[i] / np.linalg.norm(r_exact[:,i])) # compute relative error
    return t, rel_error, np.nanmax(delta)

# load analytical positions
analytical_dt_1 = load_file_mat("Analytical_dt_1.dat")
analytical_dt_01 = load_file_mat("Analytical_dt_0.1.dat")
analytical_dt_001 = load_file_mat("Analytical_dt_0.01.dat")
analytical_dt_0001 = load_file_mat("Analytical_dt_0.001.dat")
analytical_dt_00001 = load_file_mat("Analytical_dt_0.0001.dat")

# compute relative error and delta for RK4
t1, rel_error_dt1, delta1 = relative_error(pos_dt_1, analytical_dt_1)
t2, rel_error_dt01, delta2 = relative_error(pos_dt_01, analytical_dt_01)
t3, rel_error_dt001, delta3 = relative_error(pos_dt_001, analytical_dt_001)
t4, rel_error_dt0001, delta4 = relative_error(pos_dt_0001, analytical_dt_0001)
t5, rel_error_dt00001, delta5 = relative_error(pos_dt_00001, analytical_dt_00001)

# plot relative error for RK4
plt.plot(t1, rel_error_dt1, label=r"h=1")
plt.plot(t2, rel_error_dt01, label=r"h=10$^{-1}$")
plt.plot(t3, rel_error_dt001, label=r"h=10$^{-2}$")
plt.plot(t4, rel_error_dt0001, label=r"h=10$^{-3}$")
plt.plot(t5, rel_error_dt00001, label=r"h=10$^{-4}$")
plt.xlabel(r"Time [$\mu$s]")
plt.ylabel(r"Relative error")
plt.title("Runge Kutta 4")
plt.legend()
plt.yscale("log")
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("relative_error_RK4.pdf")
# plt.show()


# compute relative error and delta for forward Euler
t1, rel_error_dt1_Euler, delta1_E = relative_error(pos_dt_1_Euler, analytical_dt_1)
t2, rel_error_dt01_Euler, delta2_E = relative_error(pos_dt_01_Euler, analytical_dt_01)
t3, rel_error_dt001_Euler, delta3_E = relative_error(pos_dt_001_Euler, analytical_dt_001)
t4, rel_error_dt0001_Euler, delta4_E = relative_error(pos_dt_0001_Euler, analytical_dt_0001)
t5, rel_error_dt00001_Euler, delta5_E = relative_error(pos_dt_00001_Euler, analytical_dt_00001)

# plot relative error forward Euler
plt.plot(t1, rel_error_dt1_Euler, label=r"h=1")
plt.plot(t2, rel_error_dt01_Euler, label=r"h=10$^{-1}$")
plt.plot(t3, rel_error_dt001_Euler, label=r"h=10$^{-2}$")
plt.plot(t4, rel_error_dt0001_Euler, label=r"h=10$^{-3}$")
plt.plot(t5, rel_error_dt00001_Euler, label=r"h=10$^{-4}$")
plt.xlabel(r"Time [$\mu$s]")
plt.ylabel(r"Relative error")
plt.title("Forward Euler")
plt.legend()
plt.yscale("log")
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("relative_error_Euler.pdf")
# plt.show()

dt = np.array([1, 0.1, 0.01, 0.001, 0.0001]) # time step list

# compute error convergence rate for RK4
a = np.log10(delta2/delta1)/np.log10(dt[1]/dt[0])
b = np.log10(delta3/delta2)/np.log10(dt[2]/dt[1])
c = np.log10(delta4/delta3)/np.log10(dt[3]/dt[2])
d = np.log10(delta5/delta4)/np.log10(dt[4]/dt[3])
conv_rate_RK4 = (a+b+c+d)/4
print("Error convergence rate RK4: {:.4f}".format(conv_rate_RK4))

# compute error convergence rate for forward Euler
a_E = np.log10(delta2_E/delta1_E)/np.log10(dt[1]/dt[0])
b_E = np.log10(delta3_E/delta2_E)/np.log10(dt[2]/dt[1])
c_E = np.log10(delta4_E/delta3_E)/np.log10(dt[3]/dt[2])
d_E = np.log10(delta5_E/delta4_E)/np.log10(dt[4]/dt[3])
conv_rate_Euler = (a_E+b_E+c_E+d_E)/4
print("Error convergence rate forward Euler: {:.4f}".format(conv_rate_Euler))



"""
Look at the number of particles remaining inside the trap for different frequencies
"""

omega_v = np.arange(0.2, 2.5, 0.01) # frequency list

# create list with number of particles inside the trap
# and plot as a function of frequency
n_in = load_file_mat("particles_in_trap.dat")
plot(omega_v, n_in[0]/100, "f=0.1", r"Fraction of particles inside trap after 500 $\mu$s", r"$\omega_v$ [MHz]", r"n$_{in}$/n$_{tot}$")
plot(omega_v, n_in[1]/100, "f=0.4", r"Fraction of particles inside trap after 500 $\mu$s", r"$\omega_v$ [MHz]", r"n$_{in}$/n$_{tot}$")
plot(omega_v, n_in[2]/100, "f=0.7", r"Fraction of particles inside trap after 500 $\mu$s", r"$\omega_v$ [MHz]", r"n$_{in}$/n$_{tot}$")
plt.legend()
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("np_omega_v_all_freqs_fraction.pdf")
# plt.show()

# compute and output resonance frequency
print("Resonance frequency: {:.3f} MHz".format(omega_v[np.argmin(n_in[0])]))

# frequency list with frequencies close to one of the resonance frequencies
omega_v_new = np.arange(0.4, 0.5, 0.001)

# # particle paramters and Penning trap properties
q = 1
m = 40
B0 = 9.65*10
V0 = 0.0025*9.65*10**7
d = 0.05*10**4
omega_z = np.sqrt(2*q*V0/(m*d**2))
omega_0 = q*B0/m
omega_p = (omega_0 + np.sqrt(omega_0**2 - 2*omega_z**2))/2
omega_m = (omega_0 - np.sqrt(omega_0**2 - 2*omega_z**2))/2
print("omega_z: {:.4f} MHz".format(omega_z))
print("omega_plus: {:.4f} MHz".format(omega_p))
print("omega_minus: {:.4f} MHz".format(omega_m))


# create list with number of particles inside the trap
# and plot as a function of frequency for simulation without particle interactions
n_in_zoom_noints = load_file_mat("particles_in_trap_zoom_noints.dat")
plot(omega_v_new, n_in_zoom_noints[0]/100, "f=0.1", r"No interactions", r"$\omega_v$ [MHz]", r"n$_{in}$/n$_{tot}$")
plot(omega_v_new, n_in_zoom_noints[1]/100, "f=0.4", r"No interactions", r"$\omega_v$ [MHz]", r"n$_{in}$/n$_{tot}$")
plot(omega_v_new, n_in_zoom_noints[2]/100, "f=0.7", r"No interactions", r"$\omega_v$ [MHz]", r"n$_{in}$/n$_{tot}$")
plt.legend()
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("np_inside_trap_zoom_resfreq_noints_fraction.pdf")
# plt.show()

# create list with number of particles inside the trap
# and plot as a function of frequency for simulation with particle interactions
n_in_zoom_ints = load_file_mat("particles_in_trap_zoom_ints.dat")
plot(omega_v_new, n_in_zoom_ints[0]/100, "f=0.1", r"With interactions", r"$\omega_v$", r"n$_{in}$/n$_{tot}$")
plot(omega_v_new, n_in_zoom_ints[1]/100, "f=0.4", r"With interactions", r"$\omega_v$", r"n$_{in}$/n$_{tot}$")
plot(omega_v_new, n_in_zoom_ints[2]/100, "f=0.7", r"With interactions", r"$\omega_v$", r"n$_{in}$/n$_{tot}$")
plt.legend()
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("np_inside_trap_zoom_resfreq_ints_fraction.pdf")
# plt.show()
