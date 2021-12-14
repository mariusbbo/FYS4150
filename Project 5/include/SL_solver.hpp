#ifndef __SL_solver_hpp__
#define __SL_solver_hpp__

#include "armadillo"
#include "complex"

// Class for solving the dimensionless Schrödinger equation with the Crank-Nicolson scheme
// The constructor takes the step size h as input and initializes the following:
//  - M_: Number of points along the x- and y-axis
//  - A_, B_: Matrices used in the matrix equation for updating the wave function at each timestep
//  - u_: Matrix for storing the current wave function at each timestep
//  - V_: Matrix for the potential which corresponds to the wall with one or more slits
//  - x_, y_: Vectors with x- and y-points used to create the initial state of the wave function

// The class contains the following functions:
// SL_solver: Constructor, defines the parameters listed above.
//  Input: - Step size h

// set_potential: Creates potential V corresponding to wall at the
// center of the box. The wall can have zero, one, two or three slits.
//  Input: - dx: Thickness of wall
//         - centre_x: Position where the wall is centered
//         - dy: Length of wall piece separating the slits
//         - opening: Length of slits
//         - v0: Value of potential at position of the wall. For no wall, v0 must be set to zero
//         - n_slits: Number of slits, 0, 1, 2 or 3
//  No output, only creates potential for use inside the class

// vector_index: Converts position of matrix into position of vector.
//  Input: - i,j: Indices of position in matrix
//  Output: - k: Index of vector corresponding to the indices i,j of the matrix

// set_initial_state: Creates intial state of the wave function by computing values from
// a 2D Gaussian wave packet.
//  Input: - xc, yc: Centre of Gaussian wave packet in x- and y-direction
//        - sigma_x, sigma_y: Width of the Gaussian wave packet in x- and y-direction
//        - px, py: Momentum of the wave packet in the x- and y-direction
//  No output

// A_and_B_setup: Set up the matrices A and B that will be used to update the
// wave function at each timestep.
//  Input: - dt: Timestep
//  No output, only creates A and B for use inside the class

// solve_matrix_equation: First computes Bu^n=b, then solve the matrix equation
// Au^(n+1) = b by using the armadillo sparse matrix solver, spsolve which uses
// LU decomposition.
// No input or output, only updates the wave function matrix.

class SL_solver
{
  public:
    double h_; // Step size
    int M_; // Number of points along the x- and y-axis
    arma::sp_cx_mat A_; // Matrix A used for solving the Schrödinger equation as matrix equation
    arma::sp_cx_mat B_; // Matrix B used for solving the Schrödinger equation as matrix equation
    arma::cx_vec u_; // Matrix for initial wave function
    arma::cx_mat V_; // Matrix for potential wall
    arma::vec x_; // Vector for x points
    arma::vec y_; // Vector for y points

    // Constructor
    SL_solver(double h);

    // Create potential wall with slits
    void set_potential(double dx, double centre_x, double dy, double opening, double v0, int n_slits);

    // Convert matrix indices i,j into vector index k
    int vector_index(int i, int j);

    // Set initial state of the wave function
    void set_initial_state(double xc, double sigma_x, double px, double yc, double sigma_y, double py);

    // Create matrices A and B
    void A_and_B_setup(double dt);

    // Solve the Schrödinger equation as a matrix equation using A and B
    void solve_matrix_equation();

};

#endif
