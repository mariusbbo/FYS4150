# Project 4

Simulations of the Schrödinger equation to study the double-slit experiment. We use the Crank-Nicolson approach and take advantage of Armadillo's built-in solver which uses LU decomposition to solve the Schrödinger equation as a matrix equation.

Files:
- include/SL_solver.hpp
- src/SL_solver.cpp
- include/set_params.hpp
- src/set_params.cpp

SL_solver.hpp: Contains Class SL_solver which declares all functions used to solve the Schrödinger equation by the Crank-Nicolson approach. The functions are listed below:
-  set_potential: Creates array with high values at the positions of the wall with slits.
-  vector_index: Takes matrix indices i,j and find the corresponding vector index k.
-  set_initial_state: Creates initial state of the wave function to a quantum mechanical Gaussian wave packet.
-  A_and_B_setup: Set up the matrices A and B for solving the Schrödinger equation as a matrix equation in the Crank-Nicolson approach.
-  solve_matrix_equation: Computes matrix product Bu^n=b and solves the matrix equation Au^n+1=b. Uses the built-in Armadillo solver spsolve for solving sparse matrices by LU decomposition.

SL_solver.cpp: Take in the timestep and declares the size of the simulation region and the matrices for the wave function, the potential, A and B. Contains the functions of SL_solver.hpp.
               
set_params.hpp: Declare the function set_params which creates a dictionary containing the step size, timestep, total simulation time and parameters for the initial state of the wave function and the double-slit wall.

set_params.cpp: Contains the function set_params which sets the following parameters:
-  h        : Step size (delta x = delta y)
-  dt       : Timestep
-  T        : Total simulation time.
-  xc       : x-coordinate of centre of initial Gaussian wave packet.
-  yc       : y-coordinate of centre of initial Gaussian wave packet.
-  sigma_x  : Width of initial Gaussian wave packet in x-direction.
-  sigma_y  : Width of initial Gaussian wave packet in y-direction.
-  px       : Momentum of wave packet in x-direction.
-  py       : Momentum of wave packet in y-direction.
-  v0       : Value of potential for wall with slits. For no wall, this value must be set to zero.
-  dx       : Thickness of double-slit wall.
-  centre_x : x-position of centre of double-slit wall.
-  dy       : Length of wall piece between slits.
-  opening  : Length of slits.
-  n_slits  : Number of slits. Can be set to 0, 1, 2 or 3.




main.cpp
----
C++ code that runs creates the dictionary from the set_params function. The parameters are used to set up the initial wave function, the potential wall and the matrices A and B. Then the simulation is run by solving the Schrödinger equation in a loop over a total simulation time T for n timesteps. The real and imaginary part of the wave function on vector form are stored in matrices at each timestep. After the simulation these matrices are saved to file.

Build: g++ main.cpp src/SL_solver.cpp src/set_params.cpp -I include/ -o main.exe -larmadillo

Run: ./main.exe

plot.py
----
Python script that create all the plots for the problems in Project 5. First it converts the matrix with the wave function on vector form for each timestep into a cube with the wave function on matrix form for each timestep.

Run command: python3 plot.py

