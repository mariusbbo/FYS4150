# Project 2

Solving an eigenvalue problem for the case with a one-dimensional buckling beam using Jacobi's rotational method.

Project2.cpp
----
C++ code that does all tasks in Project 2. Contains a function creating a tridiagonal matrix called 
create_tridiagonal. Another function, called max_offdiag_symmetric, computes the maximum off-diagonal 
element of this matrix. Then the matrix and this value are used in a function called jacobi_rotate which 
performs Jacobi's rotation method. A function called jacobi_eigensolver runs jacobi_rotate until convergence is reached. 
Other functions were added to sole the problems in the project:

- analytical_eigvals_and_eigvecs: Computes eigenvalues and eigenvectors analytically, used to compare analytical results 
with numerical results.

- compare_eigvals_and_eigvecs: Compare eigenvalues and eigenvectors found with armadillo with the eigenvalues 
and eigenvectors found analytically.

- test_max_off_diag: Tests the function max_offdiag_symmetric by sending in a test matrix.

- estimate_iterations: Runs Jacobi's rotation method for different values of N to find the relation between 
the N and the number of iterations. Writes N and number of iterations to file called Nvals_iterations.txt. 
Can also run Jacobi's rotation method with dense matrix and write to the file Nvals_iterations_random.txt, 
but this has to be specified manually.

- eigvecs_diffeq: Computes eigenvalues and eigenvectors for number of steps n=10 (N=9) or n=100 (N=99) and write the 
eigenvectors corresponding to the three first eigenvalues to file. To use n=10 (N=9) or n=100 (N=99) has to be 
specified manually.

All problems in the project can be checked by running the functions mentioned above in the main function 
at the bottom of the code file. 

Build command: g++ Project2.cpp -o Project.exe -larmadillo

Run command: ./Project2.exe


Problem6.py
----
Python script that reads the file Nvals_iterations.txt, or Nvals_iterations_random if specified, 
and plot number of iterations as function of N.

Run command: python3 Problem6.py


Problem7.py
----
Python script that reads the files numerical_eigenvectors_n10.txt and analytical_eigenvectors_n10.txt if n=10.
n10 should be changed to n100 if n=100. Plot subplots comparing eigenvectors found by Jacobi's rotation method
with analytical eigenvectors for n=10 and n=100.

Run command: python3 Problem7.py
