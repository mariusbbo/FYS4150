# Project 4

Ising model simulations with the Markov chain Monte Carlo method. 

Files:
- include/MCMC.hpp
- src/MCMC.cpp
- include/IsingModel.hpp
- src/IsingModel.cpp

IsingModel.hpp: Contains class IsingModel which declares all function for the Ising model. The functions are listed below:
-  initialize_energy: Computes the energy of the initial state.
-  initialize_magnetization: Computes the magnetization of the initial state.
-  create_index_list: Creates list with indices used for determining probability ratios and change in energy. Also accounts for periodic boundary conditions.
-  get_index: Finds the index of the index list created with create_index_list which determines the probability ratio and change in energy of a particular spin flip.
-  probability_ratios: Create list with probability ratios.
-  delta_E: Create list with the possible changes in energy.
-  find_histogram_pos: Finds position of the current energy value in the binned list with energy values when creating the histograms.
-  analytical: Compute partition function, expectation values of the average energy per spin, average magnetization per spin, heat capacity and magnetic susceptibility          analytically. 
                
IsingModel.cpp: Take in lattice size L and a bool which determines if the initial state should be random or ordered. Initializes the lattice with spin states +1 and -1. The initialized lattice can be both random and ordered. Contains all functions in IsingModel.hpp.

MCMC.hpp: Contains class MCMC which declares the function run_MCMC for running the Markov chain Monte Carlo method for the Ising model.
-  run_MCMC: 
    - Contains loops over temperature, number of Monte Carlo cycles and number of attempted spin flips.
    - Create instance of IsingModel class. 
    - Runs the Markov chain Monte Carlo method for all temperatures, number of Monte Carlo cycles and number of attempted spin flips.
    - Compute and save analytical results to file.
    - Create and save histogram for energy values.
    - Compute and save average energy per spin, average magnetization per spin, heat capacity and magnetic susceptibilty to file.
                 
MCMC.cpp: Initialize number of temperature values, number of Monte Carlo cycles, lattice size, number of Monte Carlo cycles that will be ignored due to burn-in, bool which determine if initial state should be random or ordered. Contains the function run_MCMC.

main.cpp
----
C++ code that runs the run_MCMC function from the MCMC class. Run timing tests with and without parallelization. Set values for the following variables:
      - Temperature range
      - Number of temperature values
      - Number of Monte Carlo cycles
      - Lattice size
      - Number of Monte Carlo cycles that will be ignored due to burn-in
      - Bool ordered which determines if initial state is random or ordered
    
Build: 	g++ main.cpp src/IsingModel.cpp src/MCMC.cpp -I include/ -o main.exe -larmadillo -fopenmp

Run: ./main.exe

plot.py
----
Python script that create all the plots for the problems in Project 4. Also does linear regression to find the critical temperature.

Run command: python3 plot.py

