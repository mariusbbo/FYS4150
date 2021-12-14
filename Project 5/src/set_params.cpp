#include "set_params.hpp"
#include "map"


// Function creating a dictionary with the input parameters for solving
// the Schrödinger equation
// No input.
//  Output: - params: Dictionary with values of parameters for simulation
std::map<std::string, double> set_params()
{
  std::map<std::string, double> params;
  params["h"] = 0.005; // Step size (delta x = delta y = h)
  params["dt"] = 2.5e-5; // Timestep
  params["T"] = 0.002; // Total simulation time
  params["xc"] = 0.25; // x-coordinate of center of initial Gaussian wave packet
  params["sigma_x"] = 0.05; // Width of initial Gaussian wave packet in x-direction
  params["px"] = 200; // Momentum of wave packet in x-direction
  params["yc"] = 0.5; // y-coordinate of center of initial Gaussian wave packet
  params["sigma_y"] = 0.2; // Width of initial Gaussian wave packet in y-direction
  params["py"] = 0.; // Momentum of wave packet in y-direction
  params["v0"] = 1e10; // Value of potential for wall with slits. For no wall, this value must be set to zero
  params["dx"] = 0.02; // Thickness of wall
  params["centre_x"] = 0.5; // x-position of centre of wall
  params["dy"] = 0.05; // Length of wall piece between slits
  params["opening"] = 0.05; // Length of slits
  params["n_slits"] = 2.; // Number of slits. Can be set to 0, 1, 2 or 3
  return params;
}
