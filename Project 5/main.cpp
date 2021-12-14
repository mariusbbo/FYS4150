#include "SL_solver.hpp"
#include "set_params.hpp"
#include "armadillo"
#include "complex"


int main()
{
  // Dictionary with input parameters
  std::map<std::string, double> params = set_params();

  // Define instance of the SL_solver class
  SL_solver SL = SL_solver(params["h"]);

  // Set potential corresponding to wall with one or more slits
  SL.set_potential(params["dx"], params["centre_x"], params["dy"], params["opening"], params["v0"], int(params["n_slits"]));

  // Creates initial state of the wave function
  SL.set_initial_state(params["xc"], params["sigma_x"], params["px"], params["yc"], params["sigma_y"], params["py"]);

  // Creates the matrices A and B used for solving the Schrödinger equation as a matrix equation
  SL.A_and_B_setup(params["dt"]);

  int M = 1.0/params["h"]; // Number of points along the x- and y-axis

  int n = params["T"]/params["dt"] + 1; // Number of timesteps

  // Matrices for storing vectors with the real and imaginary part of the wave function
  arma::mat u_real(n, (M-2)*(M-2), arma::fill::zeros);
  arma::mat u_imag(n, (M-2)*(M-2), arma::fill::zeros);

  std::cout << "Number of timesteps: " << n << "\n";

  // Solve the Schrödinger equation and store the values of the wave function at each timestep
  for (int i = 0; i < n; i++)
  {
    std::cout << "\rTimestep: " << i+1 << std::flush;
    SL.solve_matrix_equation();
    for (int j = 0; j < (M-2)*(M-2); j++)
    {
      u_real(i,j) = real(SL.u_(j));
      u_imag(i,j) = imag(SL.u_(j));
    }
  }

  // Save the matrices containing the vectors for the real and imaginary parts of the wave function to files
  // u_real.save("u_real_doubleslit.dat");
  // u_imag.save("u_imag_doubleslit.dat");

  return 0;
}
