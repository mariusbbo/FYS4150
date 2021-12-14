#ifndef __set_params_hpp__
#define __set_params_hpp__

#include "map"

// Function setting values for input parameters and creating dictionary with these
// input parameters for simulation of the Schrödinger equation.
// No input.
//  Output: - Dictionary containing input parameters
std::map<std::string, double> set_params();

#endif
