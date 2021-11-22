import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa
import scipy.stats
plt.rcParams.update({'font.size': 16})

# function loading data and creating array
def load_file_mat(filename):
    object = pa.mat()
    object.load("{}".format(filename))
    return np.array(object)

# function for plotting
def plot(x, y, label, title, xlabel, ylabel, color):
    plt.plot(x, y, "{}".format(color), label="{}".format(label))
    plt.title("{}".format(title))
    plt.xlabel("{}".format(xlabel))
    plt.ylabel("{}".format(ylabel))


"""Compare numerical and analytical results for 2x2 lattice"""

# function which output values for epsilon, m, Cv and chi
def get_values(filename):
    num = load_file_mat(filename)
    eps = num[0][0]
    m = num[0][1]
    Cv = num[0][2]
    chi = num[0][3]
    return eps, m, Cv, chi

# defin lists for storing values
eps_num = np.zeros(5)
m_num = np.zeros(5)
Cv_num = np.zeros(5)
chi_num = np.zeros(5)

# list of files with numerical results for different number of MC cycles
filename_list = np.array(["numerical_n100.dat", "numerical_n1e3.dat", "numerical_n1e4.dat", "numerical_n1e5.dat", "numerical_n1e6.dat"])

# store values for different number of MC cycles in lists
for i in range(len(filename_list)):
    eps_num[i] = get_values(filename_list[i])[0]
    m_num[i] = get_values(filename_list[i])[1]
    Cv_num[i] = get_values(filename_list[i])[2]
    chi_num[i] = get_values(filename_list[i])[3]

# load and create arrays with analytical values
analytical = load_file_mat("analytical_results.dat")
eps_a = np.ones(5)*analytical[0][0]
m_a = np.ones(5)*analytical[0][1]
Cv_a = np.ones(5)*analytical[0][2]
chi_a = np.ones(5)*analytical[0][3]

# compute relative error for different parameters
rel_error_eps = abs((eps_num-eps_a)/eps_a)
rel_error_m = abs((m_num-m_a)/m_a)
rel_error_Cv = abs((Cv_num-Cv_a)/Cv_a)
rel_error_chi = abs((chi_num-chi_a)/chi_a)

n_list = np.array([100, 1e3, 1e4, 1e5, 1e6]) # list with number of MC cycles for plotting

# plot <epsilon>, <|m|>, C_v and chi as function of number of MC cycles

plot(n_list, eps_a, r"Analytical", r"Average energy per spin $\langle\epsilon\rangle$", r"Number of MC cycles", r"$\langle\epsilon\rangle$ [J]", "C0")
plot(n_list, eps_num, r"Numerical", r"Average energy per spin $\langle\epsilon\rangle$", r"Number of MC cycles", r"$\langle\epsilon\rangle [J]$ ", "C1")
plt.xscale("log")
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.legend()
plt.savefig("eps_an_num.pdf")
# plt.show()

plot(n_list, m_a, r"Analytical", r"Average magnetization per spin $\langle|m|\rangle$", r"Number of MC cycles", r"$\langle|m|\rangle$ ", "C0")
plot(n_list, m_num, r"Numerical", r"Average magnetization per spin $\langle|m|\rangle$", r"Number of MC cycles", r"$\langle|m|\rangle$ ", "C1")
plt.xscale("log")
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.legend()
plt.savefig("m_an_num.pdf")
# plt.show()

plot(n_list, Cv_a, r"Analytical", r"Heat capacity $C_v$", r"Number of MC cycles", r"$C_v$ [$k_B$]", "C0")
plot(n_list, Cv_num, r"Numerical", r"Heat capacity $C_v$", r"Number of MC cycles", r"$C_v$ [$k_B$]", "C1")
plt.xscale("log")
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.legend()
plt.savefig("Cv_an_num.pdf")
# plt.show()

plot(n_list, chi_a, r"Analytical", r"Magnetic susceptibility $\chi$", r"Number of MC cycles", r"$\chi$ [1/J]", "C0")
plot(n_list, chi_num, r"Numerical", r"Magnetic susceptibility $\chi$", r"Number of MC cycles", r"$\chi$ [1/J]", "C1")
plt.xscale("log")
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.legend()
plt.savefig("chi_an_num.pdf")
# plt.show()

# plot relative error for epsilon, m, Cv and chi as function of number of MC cycles

plot(n_list, rel_error_eps/1e-3, r"$\langle\epsilon\rangle$", r"Relative error for different parameters", r"Number of MC cycles", r"Relative error $\times 10^{-3}$ ", "C0")
plot(n_list, rel_error_m/1e-3, r"$\langle|m|\rangle$", r"Relative error for different parameters", r"Number of MC cycles", r"Relative error $\times 10^{-3}$", "C1")
plt.xscale("log")
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.legend()
plt.savefig("rel_error_eps_m.pdf")
# plt.show()

plot(n_list, rel_error_Cv, r"$C_v$", r"Relative error for different parameters", r"Number of MC cycles", r"Relative error", "C2")
plot(n_list, rel_error_chi, r"$\chi$", r"Relative error for different parameters", r"Number of MC cycles", r"Relative error", "C3")
plt.xscale("log")
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.legend()
plt.savefig("rel_error_cv_chi.pdf")
# plt.show()


"""Study burn-in time"""

# load files with epsilon and m for different temperatures
epsilon = load_file_mat("epsilon.dat")
eps_1 = epsilon[0]
eps_24 = epsilon[1]

m = load_file_mat("m.dat")
m_1 = m[0]
m_24 = m[1]

epsilon_o = load_file_mat("epsilon_ordered.dat")
eps_1_ordered = epsilon_o[0]
eps_24_ordered = epsilon_o[1]

m_o = load_file_mat("m_ordered.dat")
m_1_ordered = m_o[0]
m_24_ordered = m_o[1]

n = np.linspace(0, 10**6, 10**6) # list with number of MC cycles for plotting

# plot epsilon and m as function of number of MC cycles

plot(n, eps_1, r"T=1 J/$k_B$", "Random initial states", "Number of MC cycles", r"$\langle\epsilon\rangle$ [J]", "C0")
plot(n, eps_24, r"T=2.4 J/$k_B$", "Random initial states", "Number of MC cycles", r"$\langle\epsilon\rangle$ [J]", "C1")
plt.xscale("log")
figure = plt.gcf()
figure.set_size_inches(9, 6)
plt.legend()
plt.savefig("eps_random.pdf")
# plt.show()

plot(n, eps_1_ordered, r"T=1 J/$k_B$", "Ordered initial states", "Number of MC cycles", r"$\langle\epsilon\rangle$ [J]", "C0")
plot(n, eps_24_ordered, r"T=2.4 J/$k_B$", "Ordered initial states", "Number of MC cycles", r"$\langle\epsilon\rangle$ [J]", "C1")
plt.xscale("log")
figure = plt.gcf()
figure.set_size_inches(9, 6)
plt.legend()
plt.savefig("eps_ordered.pdf")
# plt.show()

plot(n, m_1, r"T=1 J/$k_B$", r"Random initial states", "Number of MC cycles", r"$\langle|m|\rangle$", "C0")
plot(n, m_24, r"T=2.4 J/$k_B$", r"Random initial states", "Number of MC cycles", r"$\langle|m|\rangle$", "C1")
plt.xscale("log")
figure = plt.gcf()
figure.set_size_inches(9, 6)
plt.legend()
plt.savefig("m_random.pdf")
# plt.show()

plot(n, m_1_ordered, r"T=1 J/$k_B$", r"Ordered initial states", "Number of MC cycles", r"$\langle|m|\rangle$", "C0")
plot(n, m_24_ordered, r"T=2.4 J/$k_B$", r"Ordered initial states", "Number of MC cycles", r"$\langle|m|\rangle$", "C1")
plt.xscale("log")
figure = plt.gcf()
figure.set_size_inches(9, 6)
plt.legend()
plt.savefig("m_ordered.pdf")
# plt.show()


"""Create histograms"""

# load and create array with values for histogram
histogram = load_file_mat("E_histogram.dat")
eps_hist_T1 = histogram[0] # values for temperature T=1
eps_hist_T24 = histogram[1] # values for temperature T=1

eps = np.linspace(-4, 4, len(eps_hist_T24)) # epsilon list for plotting

# plot normalised histograms for epsilon for T=1 and T=2.4

plt.bar(eps, eps_hist_T1/400/sum(eps_hist_T1)/10**-4, width=8/400, align='edge')
plt.title("Normalised p($\epsilon$;T), T=1 [J/$k_B$]")
plt.xlabel("$\epsilon$ [J]")
plt.ylabel(r"p($\epsilon$;T=1) $\times 10^{-4}$")
plt.xlim(-2.1, -1.8)
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.savefig("pdf_T1.pdf")
# plt.show()

plt.bar(eps, eps_hist_T24/400/sum(eps_hist_T24)/10**-4, width=8/400, align='edge')
plt.title("Normalised p($\epsilon$;T), T=2.4 [J/$k_B$]")
plt.xlabel("$\epsilon$ [J]")
plt.ylabel(r"p($\epsilon$;T=2.4) $\times 10^{-4}$")
plt.xlim(-2.5, 0)
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.savefig("pdf_T24.pdf")
# plt.show()


"""Timing tests"""

# load and create files for comparing simulation time with and without parallelization
params_para = load_file_mat("n_time_list_para.dat")
n = params_para[0]
time_para = params_para[1]

params_nopara = load_file_mat("n_time_list_no_para.dat")
time_nopara = params_nopara[1]

# plot simulation time as function of number of cycles

#def plot(x, y, label, title, xlabel, ylabel, color):

plot(n, time_para, "With parallelization", "Timing tests with L = 20", "Number of MC cycles", "Time [seconds]", "C0")
plot(n, time_nopara, "No parallelization", "Timing tests with L = 20", "Number of MC cycles", "Time [seconds]", "C1")
plt.xscale("log")
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.legend()
plt.grid(alpha=0.3)
plt.savefig("timing_tests.pdf")
# plt.show()


"""Investigating phase transitions"""

# load files and create arrays for <epsilon>, <|m|>, C_v and chi

mat_L40 = load_file_mat("params_L40.dat")
eps_L40 = mat_L40[:,0]
m_L40 = mat_L40[:,1]
Cv_L40 = mat_L40[:,2]
chi_L40 = mat_L40[:,3]

mat_L60 = load_file_mat("params_L60.dat")
eps_L60 = mat_L60[:,0]
m_L60 = mat_L60[:,1]
Cv_L60 = mat_L60[:,2]
chi_L60 = mat_L60[:,3]

mat_L80 = load_file_mat("params_L80.dat")
eps_L80 = mat_L80[:,0]
m_L80 = mat_L80[:,1]
Cv_L80 = mat_L80[:,2]
chi_L80 = mat_L80[:,3]

mat_L100 = load_file_mat("params_L100.dat")
eps_L100 = mat_L100[:,0]
m_L100 = mat_L100[:,1]
Cv_L100 = mat_L100[:,2]
chi_L100 = mat_L100[:,3]

T = np.linspace(2.1, 2.4, 100) # temperature list for plotting
T_L80 = np.linspace(2.1, 2.4, 50) # temperature list for plotting
L_list = np.array([40, 60, 80, 100]) # lattice size list

# critical temperature Tc list for different lattice sizes L
Tc_list = np.array([T[np.argmax(Cv_L40)], T[np.argmax(Cv_L60)], T[np.argmax(Cv_L80)], T[np.argmax(Cv_L100)]])

# perform linear regression to determine Tc(L=inf)
T_critical = scipy.stats.linregress(L_list, Tc_list)[1] # critical temperature

# output Tc(L=inf)
# print("Critical temperature values T_c = {:.5f}".format(T_critical))

# plot <epsilon>, <|m|>, Cv and chi as function of temperature
# for different lattice sizes

plot(T, eps_L40, "L = 40", r"Average energy per spin $\langle\epsilon\rangle", r"Temperature [J/$k_B$]", r"$\langle\epsilon\rangle$ [J]", "C0")
plot(T, eps_L60, "L = 60", r"Average energy per spin $\langle\epsilon\rangle$", r"Temperature [J/$k_B$]", r"$\langle\epsilon\rangle$ [J]", "C1")
plot(T, eps_L80, "L = 80", r"Average energy per spin $\langle\epsilon\rangle$", r"Temperature [J/$k_B$]", r"$\langle\epsilon\rangle$ [J]", "C2")
plot(T, eps_L100, "L = 100", r"Average energy per spin $\langle\epsilon\rangle$", r"Temperature [J/$k_B$]", r"$\langle\epsilon\rangle$ [J]", "C3")
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.legend()
plt.savefig("eps_T.pdf")
# plt.show()

plot(T, m_L40, "L = 40", r"Average magnetization per spin $\langle|m|\rangle$", r"Temperature [J/$k_B$]", r"$\langle|m|\rangle$", "C0")
plot(T, m_L60, "L = 60", r"Average energy per spin $\langle|m|\rangle$", r"Temperature [J/$k_B$]", r"$\langle|m|\rangle$", "C1")
plot(T, m_L80, "L = 80", r"Average energy per spin $\langle|m|\rangle$", r"Temperature [J/$k_B$]", r"$\langle|m|\rangle$", "C2")
plot(T, m_L100, "L = 100", r"Average energy per spin $\langle|m|\rangle$", r"Temperature [J/$k_B$]", r"$\langle|m|\rangle$", "C3")
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.legend()
plt.savefig("m_T.pdf")
# plt.show()

plot(T, Cv_L40, "L = 40", r"Heat capacity $C_v$", r"Temperature [J/$k_B$]", r"$C_v$ [$k_B$]", "C0")
plot(T, Cv_L60, "L = 60", r"Heat capacity $C_v$", r"Temperature [J/$k_B$]", r"$C_v$ [$k_B$]", "C1")
plot(T, Cv_L80, "L = 80", r"Heat capacity $C_v$", r"Temperature [J/$k_B$]", r"$C_v$ [$k_B$]", "C2")
plot(T, Cv_L100, "L = 100", r"Heat capacity $C_v$", r"Temperature [J/$k_B$]", r"$C_v$ [$k_B$]", "C3")
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.legend()
plt.savefig("Cv_T.pdf")
# plt.show()

plot(T, chi_L40, "L = 40", r"Magnetic susceptibility $\chi$", r"Temperature [J/$k_B$]", r"$\chi$ [1/J]", "C0")
plot(T, chi_L60, "L = 60", r"Magnetic susceptibility $\chi$", r"Temperature [J/$k_B$]", r"$\chi$ [1/J]", "C1")
plot(T, chi_L80, "L = 80", r"Magnetic susceptibility $\chi$", r"Temperature [J/$k_B$]", r"$\chi$ [1/J]", "C2")
plot(T, chi_L100, "L = 100", r"Magnetic susceptibility $\chi$", r"Temperature [J/$k_B$]", r"$\chi$ [1/J]", "C3")
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.legend()
plt.savefig("chi_T.pdf")
# plt.show()
