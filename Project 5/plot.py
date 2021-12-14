import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pyarma as pa
plt.rcParams.update({'font.size': 16})


# Function loading data and creating array
# Input: - Filename of file containing matrix
# Output: - Array with elements from input file
def load_file_mat(filename):
    object = pa.mat()
    object.load("{}".format(filename))
    return np.array(object)

# Function for plotting
# Input: - x-axis list
#        - y-axis list
#        - Label
#        - Title
#        - x-axis label
#        - y-axis label
# No output, just creates plot
def plot(x, y, label, title, xlabel, ylabel):
    plt.plot(x, y, label="{}".format(label))
    plt.title("{}".format(title))
    plt.xlabel("{}".format(xlabel))
    plt.ylabel("{}".format(ylabel))

# Function for loading files and converting input matrix with u-vectors into cube with u-matrices.
# Also computes the probability u*u
# Input: - filename_real: Filename of file containing real part of wave function u
#        - filename_imag: Filename of file containing imaginary part of wave function u
# Output: - p: Cube with probability matrices for each timestep
#         - u_cube: Cube with u-matrices for each timestep
def create_u_cube(filename_real, filename_imag):
    z = complex(0,1) # Complex number i that will be multiplied by the imaginary part of the wave function
    u_mat_real = load_file_mat(filename_real)
    u_mat_imag = load_file_mat(filename_imag)*z
    u_mat = u_mat_real + u_mat_imag

    n = len(u_mat) # Number of timesteps
    M = int(np.sqrt(len(u_mat[0]))) # Number of points along x- and y-axis (actually M-2, just called M here)

    u_cube_temp = u_mat.reshape(n,M,M) # Convert matrix with u-vector into cube with u-matrices
    u_cube = np.zeros([n,M+2,M+2], dtype=complex) # Add zeros at the boundaries
    u_cube[:,1:M+1, 1:M+1] = u_cube_temp # Cube with u-matrices including boundary conditions

    prob = np.conj(u_cube)*u_cube # Probability u*u
    p = prob.real # Get real part since prob is an imaginary cube even though the imaginary part is zero
    return p, u_cube


"""Probability deviation (Problem 7)"""

# Load files and compute probability for a box with no slits and for box with double-slit
p_nowall = create_u_cube("u_real_nowall_new.dat", "u_imag_nowall_new.dat")[0]
p_doubleslit = create_u_cube("u_real_doubleslit_new.dat", "u_imag_doubleslit_new.dat")[0]

# Arrays for storing deviations of the total probability from 1
deviation_nowall = np.zeros(len(p_nowall))
deviation_doubleslit = np.zeros(len(p_doubleslit))

# Computes the difference between the total probability and 1
for i in range(len(p_nowall)):
    deviation_nowall[i] = 1-np.sum(p_nowall[i])
    deviation_doubleslit[i] = 1-np.sum(p_doubleslit[i])

dt = 2.5e-5 # Timestep
t_points_p7 = np.arange(0, 0.008+dt, dt) # Array of time points

# Plots the deviation of the total probability from 1 as function of time for a box with no slit
plot(t_points_p7/1e-3, deviation_nowall/1e-15, "", r"Deviation of total probability from 1, no wall",
     r"Time $\times$ $10^{-3}$", r"Deviation $\times$ $10^{-15}$")
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.savefig("deviation_nowall.pdf")
plt.show()

# Plot the deviation of the total probability from 1 as function of time for a box with a double-slit
plot(t_points_p7/1e-3, deviation_doubleslit/1e-15, "", r"Deviation of total probability from 1, d   ouble-slit",
     r"Time $\times$ $10^{-3}$", r"Deviation$\times$ $10^{-15}$")
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.savefig("deviation_doubleslit.pdf")
# plt.show()


"""Colourmaps (Problem 8)"""

# Function for plotting, saving and visualising colourmap
# Use code from project description
# Input: - figure_filename: Filename of colourmap that will be saved as pdf
#        - p: Probability cube with probability grids for each timestep
#        - t_points: Array with time points
#        - t: Simulation time when colourmap is created
#        - cbar_label: Label for colourbar
#        - title: Title of colourmap that is created
# No output, just creates and save colourmap
def plot_colourmap(figure_filename, p, t_points, t, cbar_label, title):
    # Some settings
    fontsize = 16

    # Create figure
    fig = plt.figure()
    ax = plt.gca()
    figure = plt.gcf()
    figure.set_size_inches(8, 6)

    time_index = np.argwhere(t_points==t)[0,0] # Find time index for when to create colourmap

    # Create a colour scale normalization according to the max p value in the time frame when the colourmap is plotted
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(p[time_index]))

    # Create colourmap
    img = ax.imshow(p[time_index], extent=[0,1,0,1], cmap=plt.get_cmap("viridis"), norm=norm)

    # Axis labels
    plt.xlabel("x", fontsize=fontsize)
    plt.ylabel("y", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    # Add a colourbar
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label("{}".format(cbar_label), fontsize=fontsize, labelpad=15)
    cbar.ax.tick_params(labelsize=fontsize)

    # Add a text element showing the time
    time_txt = plt.text(0.95, 0.95, "t = {:.1e}".format(t), color="white",
                        horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

    ax.set_title("{}".format(title))
    plt.savefig("{}".format(figure_filename))
    # plt.show()

t_points_p8 = np.arange(0, 0.002+dt, dt) # Array with time points

# load and create probability cube and wave function cube for box with double-slit
p_doubleslit_p8 = create_u_cube("u_real_doubleslit_p8.dat", "u_imag_doubleslit_p8.dat")[0]
u_doubleslit_p8 = create_u_cube("u_real_doubleslit_p8.dat", "u_imag_doubleslit_p8.dat")[1]

# plot colourmaps of the probability for box with double-slit at different times
plot_colourmap("colourmap_t0.pdf", p_doubleslit_p8, t_points_p8, 0,
               r"$p$($x,y;t$=0) / max( $p$($x,y;t$=0) )", "Probability")
plot_colourmap("colourmap_t001.pdf", p_doubleslit_p8, t_points_p8, 0.001,
               r"$p$($x,y;t$=0.001) / max ($p$($x,y;t$=0.001) )", "Probability")
plot_colourmap("colourmap_t002.pdf", p_doubleslit_p8, t_points_p8, 0.002,
               r"$p$($x,y;t$=0.002 )/ max( $p$($x,y;t$=0.002) )", "Probability")

# plot colourmaps of the real part of the wave function for box with double-slit at different times
plot_colourmap("colourmap_t0_real.pdf", u_doubleslit_p8.real, t_points_p8, 0,
               r"Re $u$($x,y;t$=0) / max( Re $u$($x,y;t$=0) )", r"Real part of wave function")
plot_colourmap("colourmap_t001_real.pdf", u_doubleslit_p8.real, t_points_p8, 0.001,
               r"Re $u$($x,y;t$=0.001) / max( Re $u$($x,y;t$=0.001) )", r"Real part of wave function")
plot_colourmap("colourmap_t002_real.pdf", u_doubleslit_p8.real, t_points_p8, 0.002,
               r"Re $u$($x,y;t$=0.002) / max( Re $u$($x,y;t$=0.002) )", r"Real part of wave function")

# plot colourmaps of the imaginary part of the wave function for box with double-slit at different times
plot_colourmap("colourmap_t0_imag.pdf", u_doubleslit_p8.imag, t_points_p8, 0,
               r"Im $u$($x,y;t$=0) / max( Im $u$($x,y;t$=0) )", r"Imaginary part of wave function")
plot_colourmap("colourmap_t001_imag.pdf", u_doubleslit_p8.imag, t_points_p8, 0.001,
               r"Im $u$($x,y;t$=0.001) / max( Im $u$($x,y;t$=0.001) )", r"Imaginary part of wave function")
plot_colourmap("colourmap_t002_imag.pdf", u_doubleslit_p8.imag, t_points_p8, 0.002,
               r"Im $u$($x,y;t$=0.002) / max( Im $u$($x,y;t$=0.002) )", r"Imaginary part of wave function")


"""Detection probability (Problem 9)"""

h = 0.005 # step size in x- and y-direction
x = np.arange(0, 1, h) # Array with x-points
y = np.arange(0, 1, h) # Array with y-points

time_index = np.argwhere(t_points_p8==0.002)[0,0] # Find index of t=0.002 in time array
x_index = np.argmin(abs(x-0.8)) # Find index of x=0.8 in x array

# load and create probability cubes for box with a single slit and triple-slit
p_singleslit = create_u_cube("u_real_singleslit.dat", "u_imag_singleslit.dat")[0]
p_tripleslit = create_u_cube("u_real_tripleslit.dat", "u_imag_tripleslit.dat")[0]

# Find probability distribution at x=0.8 and t=0.002 for box with single slit, double-slit and triple-slit
screen_prob_doubleslit = p_doubleslit_p8[time_index,:,x_index]
screen_prob_singleslit = p_singleslit[time_index,:,x_index]
screen_prob_tripleslit = p_tripleslit[time_index,:,x_index]

# Plot the normalised probability distribution at x=0.8 and t=0.002 for box with double-slit
plot(y, screen_prob_doubleslit/sum(screen_prob_doubleslit), "", "Detection probability double-slit", "y", "p(y|x=0.8;t=0.002)")
figure = plt.gcf()
figure.set_size_inches(9, 6)
plt.savefig("screen_prob_doubleslit.pdf")
# plt.show()

# Plot the normalised probability distirbution at x=0.8 and t=0.002 for box with a single slit
plot(y, screen_prob_singleslit/sum(screen_prob_singleslit), "", "Detection probability single-slit", "y", "p(y|x=0.8;t=0.002)")
figure = plt.gcf()
figure.set_size_inches(9, 6)
plt.savefig("screen_prob_singleslit.pdf")
# plt.show()

# Plot the normalised probability distirbution at x=0.8 and t=0.002 for box with triple-slit
plot(y, screen_prob_tripleslit/sum(screen_prob_tripleslit), "", "Detection probability triple-slit", "y", "p(y|x=0.8;t=0.002)")
figure = plt.gcf()
figure.set_size_inches(9, 6)
plt.savefig("screen_prob_tripleslit.pdf")
# plt.show()


# Code for creating animations, taken from the project description
# Created function taking probability and array with time points as input to easier
# animate results from different simulations
def make_animation(p, t_points):

    # Some settings
    fontsize = 12
    t_min = t_points[0]

    # Create figure
    fig = plt.figure()
    ax = plt.gca()

    # Create a colour scale normalization according to the max p value in the first frame
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(p[0]))

    # Plot the first frame
    img = ax.imshow(p[0], extent=[0,1,0,1], cmap=plt.get_cmap("viridis"), norm=norm)

    # Axis labels
    plt.xlabel("x", fontsize=fontsize)
    plt.ylabel("y", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    # Add a colourbar
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label("p(x,y,t)", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    # Add a text element showing the time
    time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white",
                        horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

    # Function that takes care of updating the z data and other things for each frame
    def animation(i):
        # Normalize the colour scale to the current frame?
        norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(p[i]))
        img.set_norm(norm)

        # Update z data
        img.set_data(p[i])

        # Update the time label
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3e}".format(current_time))

        return img

    # Use matplotlib.animation.FuncAnimation to put it all together
    anim = FuncAnimation(fig, animation, interval=100, frames=np.arange(0, len(p), 2), repeat=False, blit=0)

    # Run the animation!
    plt.show()

# Create animations for box with single slit, double-slit, triple-slit and no slits
# Also use arrays with different time points
make_animation(np.sqrt(p_doubleslit), t_points_p7)
make_animation(p_doubleslit, t_points_p7)
make_animation(p_singleslit, t_points_p8)
make_animation(p_tripleslit, t_points_p8)
