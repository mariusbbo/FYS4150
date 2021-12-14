#include "SL_solver.hpp"
#include "armadillo"
#include "complex"

// For using complex numbers
using namespace std::complex_literals;

// Constructor
SL_solver::SL_solver(double h)
{
  h_ = h; // Step size
  M_ = 1.0/h_; // Number of points along the x- and y-axis
  x_ = arma::linspace(0, 1, M_); // Vector with x-values
  y_ = arma::linspace(0, 1, M_); // Vector with y-values
  V_ = arma::cx_mat(M_-2, M_-2, arma::fill::zeros); // Complex matrix for storing the potential
  A_ = arma::sp_cx_mat((M_-2)*(M_-2), (M_-2)*(M_-2)); // Complex sparse matrix for storing the values of A
  B_ = arma::sp_cx_mat((M_-2)*(M_-2), (M_-2)*(M_-2)); // Complex sparse matrix for storing the values of B
  u_ = arma::cx_vec((M_-2)*(M_-2), arma::fill::zeros); // Complex vector for storing the wave function at each timestep
}


// Function for setting the potential corresponding to the wall with one or more slits
//  Input: - dx: Wall thickness
//         - centre_x: x-position at which to place the wall
//         - dy: Length of wall piece between slits
//         - opening: Length of slits
//         - v0: Value of potential at position of the wall. For no wall, v0 must be set to zero
//         - n_slits: Number of slits, 0, 1, 2 or 3
//  No output
void SL_solver::set_potential(double dx, double centre_x, double dy, double opening, double v0, int n_slits)
{
  std::cout << "Number of slits: " << n_slits << "\n";
  // Creates wall with no slits if n_slits=0
  if (n_slits == 0)
  {
    int indx1 = arma::index_min(arma::abs(x_ - (centre_x - dx/2.0))); // Index of left side of the wall
    int indx2 = arma::index_min(arma::abs(x_ - (centre_x + dx/2.0))); // Index of right side of the wall
    V_.submat(indx1, 0, indx2, M_-3).fill(v0);
  }

  // Creates wall with one slit if n_slits=1
  else if (n_slits == 1)
  {
    int indx1 = arma::index_min(arma::abs(x_ - (centre_x - dx/2.0))); // Index of left side of the wall
    int indx2 = arma::index_min(arma::abs(x_ - (centre_x + dx/2.0))); // Index of right side of the wall
    int indy1 = arma::index_min(arma::abs(y_ - (0.5 - opening/2.0))); // Index of lowest y-position of slit
    int indy2 = arma::index_min(arma::abs(y_ - (0.5 + opening/2.0))); // Index of highest y-position of slit
    V_.submat(indx1, 0, indx2, indy1).fill(v0); // Create wall below the slit
    V_.submat(indx1, indy2, indx2, M_-3).fill(v0); // Create wall above the slit
  }

  // Creates wall with two slits if n_slits=2
  else if (n_slits == 2)
  {
    int indx1 = arma::index_min(arma::abs(x_ - (centre_x - dx/2.0))); // Index of left side of the wall
    int indx2 = arma::index_min(arma::abs(x_ - (centre_x + dx/2.0))); // Index of right side of the wall
    int indy1 = arma::index_min(arma::abs(y_ - (0.5 - dy/2.0 - opening))); // Index of lowest y-position of the lower slit
    int indy2 = arma::index_min(arma::abs(y_ - (0.5 - dy/2.0))); // Index of highest y-position of the lower slit
    int indy3 = arma::index_min(arma::abs(y_ - (0.5 + dy/2.0))); // Index of lowest y-position of the upper slit
    int indy4 = arma::index_min(arma::abs(y_ - (0.5 + dy/2.0 + opening))); // Index of highest y-position of the upper slit

    V_.submat(indx1, 0, indx2, indy1).fill(v0); // Create wall below the lower slit
    V_.submat(indx1, indy2, indx2, indy3).fill(v0); // Create wall between the slits
    V_.submat(indx1, indy4, indx2, M_-3).fill(v0); // Create wall above the upper slit
  }

  // Creates wall with three slits if n_slits=3
  else if (n_slits==3)
  {
    int indx1 = arma::index_min(arma::abs(x_ - (centre_x - dx/2.0))); // Index of left side of the wall
    int indx2 = arma::index_min(arma::abs(x_ - (centre_x + dx/2.0))); // Index of right side of the wall
    int indy1 = arma::index_min(arma::abs(y_ - (0.5 - 1.5*opening - dy))); // Index of the lowest y-position of the lower slit
    int indy2 = arma::index_min(arma::abs(y_ - (0.5 - opening/2.0 - dy))); // Index of the highest y-position of the lower slit
    int indy3 = arma::index_min(arma::abs(y_ - (0.5 - opening/2.0))); // Index of the lowest y-position of the middle slit
    int indy4 = arma::index_min(arma::abs(y_ - (0.5 + opening/2.0))); // Index of the highest y-position of the middle slit
    int indy5 = arma::index_min(arma::abs(y_ - (0.5 + opening/2.0 + dy))); // Index of the lowest y-position of the upper slit
    int indy6 = arma::index_min(arma::abs(y_ - (0.5 + 1.5*opening + dy))); // Index of the highest y-position of the upper slit

    V_.submat(indx1, 0, indx2, indy1).fill(v0); // Create wall below the lower slit
    V_.submat(indx1, indy2, indx2, indy3).fill(v0); // Create wall between the lower and the middle slit
    V_.submat(indx1, indy4, indx2, indy5).fill(v0); // Create wall between the middle and the upper slit
    V_.submat(indx1, indy6, indx2, M_-3).fill(v0); // Create wall above the upper slit
  }
  else
  {
    std::cout << "This code does not create more than 3 slits!" << "\n";
    std::exit(0);
  }
}


// Function for converting matrix indices i,j into corresponding vector index k
//  Input: - i,j: Matrix indices
//  Output: - k: Vector index
int SL_solver::vector_index(int i, int j)
{
  // Find index k from i and j
  return j*(M_-2) + i;
}


// Function for setting the initial state of the wave function.
// Uses a quantum mechanical Gaussian wave packet
//  Input: - xc, yc: Central x- and y-position of the wave packet
//         - sigma_x, sigma_y: Width of the wave packet in x- and y-direction
//         - px, py: Momentum in the x- and y-direction
//  No output
void SL_solver::set_initial_state(double xc, double sigma_x, double px, double yc, double sigma_y, double py)
  {
    arma::cx_mat u(M_-2, M_-2); // Define complex matrix for the wave function u

    // Loops for computing the value of the Gaussian wave packet at each point in the wave function matrix
    for (int i = 0; i < M_-2; i++)
    {
      for (int j = 0; j < M_-2; j++)
      {
        // Computes values of the Gaussian wave packet
        u(i,j) = std::exp(-(((x_(i)-xc)*(x_(i)-xc))/(2*sigma_x*sigma_x)) - (((y_(j)-yc)*(y_(j)-yc))/(2*sigma_y*sigma_y))
         + 1.0i*px*(x_(i)-xc) + 1.0i*py*(y_(j)-yc));
      }
    }
    // Find normalization constant
    arma::cx_double norm = arma::accu(arma::conj(u)%u);

    // Convert wave function matrix into a vector and normalise
    for (int i = 0; i < M_-2; i++)
    {
      for (int j = 0; j < M_-2; j++)
      {
        int k = vector_index(i,j);
        u_(k) = u(i,j)/std::sqrt(norm);
      }
    }
}


// Function for setting up the matrices A and B which is used in the matrix equation.
//  Input: - dt: Timestep
//  No output
void SL_solver::A_and_B_setup(double dt)
{
  arma::cx_double r = 1.0i*dt/(2*h_*h_); // Ratio of timestep and step size

  arma::cx_vec a((M_-2)*(M_-2), arma::fill::zeros); // Vectors for storing values of main diagonal of A
  arma::cx_vec b((M_-2)*(M_-2), arma::fill::zeros); // Vectors for storing values of main diagonal of B

  for (int i = 0; i < M_-2; i++)
  {
    for (int j = 0; j < M_-2; j++)
    {
      int k = vector_index(i,j); // Find vector index k corresponding to matrix indices i,j
      a(k) = 1.0 + 4.0*r + 1.0i*dt*V_(i,j)/2.0; // Add values to a
      b(k) = 1.0 - 4.0*r - 1.0i*dt*V_(i,j)/2.0; // Add values to b
    }
  }

  A_.diag() = a; // Set main diagonal of A to a
  A_.diag(1) += -r; // Set upper diagonal to the ratio -r
  A_.diag(-1) += -r; // Set lower diagonal to the ratio -r
  A_.diag(M_-2) += -r; // Set the upper diagonal on the x-position M-2 to the ratio -r
  A_.diag(-(M_-2)) += -r; // Set the lower diagonal on the x-position M-2 to the ratio -r

  // Set every M-2 values in the M-2 diagonals to zero
  for (int i = M_-2; i < (M_-2)*(M_-2); i += M_-2)
  {
    A_(i, i-1) = 0;
    A_(i-1, i) = 0;
  }
  // Set values of B
  B_ = A_ * -1.0;
  B_.diag() = b;
}


// Function for solving the matrix equation to update the wave function at every timestep.
// Uses armadillo's sparse matrix equation solver spsolve which uses LU decomposition.
// No input or output.
void SL_solver::solve_matrix_equation()
{
  arma::cx_vec b = B_ * u_; // Computes the vector b
  u_ = arma::spsolve(A_, b); // Update the wave function u by solving the matrix equation Au^(n+1)=b
}
