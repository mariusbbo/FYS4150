# Project 3

Penning trap simulation. 

Files:
- include/PenningTrap.hpp
- src/PenningTrap.cpp
- include/Particle.hpp
- src/Particle.cpp

Particle.hpp contains a class Particle which stores the particle parameters charge, mass, position and velocity in a list.

Particle.cpp contains the constructor of the Particle class in Particle.hpp.

PenningTrap.hpp contains a class PenningTrap which declares all functions needed to simulate a Penning trap. The functions are listed below.

  - add_particle: Adds a particle with parameters defined in the Particle class to the Penning trap.
  - external_E_field: Computes the external electric field.
  - external_B_field: Computes the external magnetic field.
  - force_particle: Computes the Coulomb force from one particle.
  - total_force_external: Computes the Lorentz force from the electric and magnetic field.
  -  total_force_particles: Computes the Coulomb force on one particle from all other particles.
  -  total_force: Adds the Coulomb force and Lorentz force computed with total_force_particle and total_force_external
  -  evolve_RK4: Evolves the system one timestep using Runge-Kutta 4.
  -  evolve_forward_Euler: Evolves the system one timestep using forward Euler.
  -  simulation: Runs evolve_RK4 or evolve_forward_Euler a specifies number of timesteps for a specified number of particles. Saves positions and velocities to file.
  -  analytical_solution: Computes analytical solution. Saves analytical positions to file.
  -  count_particles: Counts the number of particles remaining inside the trap after a specific amount of time.

PenningTrap.cpp initializes the Penning trap properties and contains the functions in the PenningTrap class.

main.cpp
----
C++ code that runs the simulations of the Penning trap. The particle parameters and the Penning trap properties are defined. 
Simulate the Penning trap for different values of the amplitude and angular frequency of the time-dependent electric potential.
Run the function count_particles and save the number of remaining particles to file.

Build: g++ main.cpp src/Particle.cpp src/PenningTrap.cpp -I include/ -o main.exe -O3 -larmadillo

Run: ./main.exe

plot.py
----
Python script that create all the plots for the problems in Project 3. It also computes the angular 
frequencies, relative errors and the error convergence rates.

Run command: python3 plot.py
