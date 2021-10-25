# Project 3

Penning trap simulation. 

Files:
- include/PenningTrap.hpp
- src/PenningTrap.cpp
- include/Particle.hpp
- src/Particle.cpp

Particle.hpp contains a class Particle which store the particle parameters charge, mass, position and velocity in a list.

Particle.cpp contains the constructor of the Particle class in Particle.hpp.

PenningTrap.hpp contains a class PenningTrap which declare all functions needed to simulate a Penning trap.

PenningTrap.cpp initializes the Penning trap properties. It is also where all the functions for the class PenningTrap 
are written. In this file the electric and magnetic field are computed and used to find the force on the particles. 
The total force on the particles are found and used in the Runge-Kutta 4 and forward Euler functions. These methods 
are simulated with a function simulate.

Build: g++ main.cpp src/Particle.cpp src/PenningTrap.cpp -I include/ -o main.exe -O3 -larmadillo

Run: ./main.exe

plot.py
----
Python script that create all the plots for the problems in Project 3. It also computes the angular 
frequencies, relative errors and the error convergence rates.

Run command: python3 plot.py
