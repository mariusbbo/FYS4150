#include "IsingModel.hpp"
#include "armadillo"

// Constructor
// Input: - lattice size L
//       - ordered, bool which determine if initial state is ordered or not
IsingModel::IsingModel(int L, bool ordered)
{
  L_ = L; // lattice size
  // initialize ordered or random state
  if (ordered==true)
  {
    S_ = arma::mat(L_, L_, arma::fill::ones);
  }
  else
  {
    arma::arma_rng::set_seed_random();
    arma::mat S = arma::randi<arma::mat>(L_, L_, arma::distr_param(0, 1)); // create lattice with 0 and 1
    S_ = 2*S - 1; // convert lattice values from 0,1 to -1,+1
  }
}

// Initialize energy
// Output: - Energy of initial state
double IsingModel::initialize_energy()
{
  // shift lattice up, down, left and right to sum over all neighbouring pairs
  arma::mat S_left = arma::shift(S_, -1, 1);
  arma::mat S_right = arma::shift(S_, 1, 1);
  arma::mat S_up = arma::shift(S_, -1, 0);
  arma::mat S_down = arma::shift(S_, 1, 0);
  arma::mat S_sum = S_%S_left + S_%S_right + S_%S_up + S_%S_down;
  return -arma::accu(S_sum);
}

// Initialize magnetization
// Output: - Magnetization of initial state
double IsingModel::initialize_magnetization()
{
  return arma::accu(S_); // sum over all spins for initial state
}

// Create list with indices which will be used to determine the change
// in energy and the probability ratio for each attempted spin flip
// The list also accounts for the periodic boundary conditions
// For a 2x2 lattice the list will look like (2, 0, 1, 2, 0)
void IsingModel::create_index_list()
{
  // defines list with two more points than the lattice size
  arma::vec i_list = arma::vec(L_+2, arma::fill::zeros);

  for (int i = -1; i < L_+1; i++)
  {
    i_list(i+1) = (L_+i) % L_;
  }
  i_list_ = i_list;
}

// Finds index from list in create_index_list which determines the change in energy and
// probability ratio
// Input: - Indices for spin that will be attempeted to be flipped, i,j
// Output: - Index of index list
int IsingModel::get_index(int i, int j)
{
  // find sum of surrounding spins
  // multiply by the value of spin to account for it being negative
  double sum_spins = (S_(i, i_list_(j)) + S_(i, i_list_(j+2)) + S_(i_list_(i), j) + S_(i_list_(i+2), j))*S_(i,j);
  return (sum_spins + 4) / 2;
}

// Create list with probability ratios
// Input: - Temperature T
// Output: - List with probability ratios
arma::vec IsingModel::probability_ratios(double T)
{
  return arma::vec({std::exp(8/T), std::exp(4/T), 1, std::exp(-4/T), std::exp(-8/T)});
}

// Create list with possible changes in energy
// Output: - List with change in energy
arma::vec IsingModel::delta_E()
{
  return arma::vec({-8, -4, 0, 4, 8});
}

// Find index of list for creating histogram of energy values
// Input: - List with energy values from -E to +E, comp_list
//        - List with ones, test
//        - Energy from current MC cycle
// Output: - Index of histogram list to which the energy will be added
int IsingModel::find_histogram_pos(arma::vec comp_list, arma::vec test, double E)
{
  arma::vec diff = arma::abs(comp_list - test*E);
  return arma::index_min(diff);
}

// Compute analytical values for epsilon, m, Cv and chi
// Input: - Temperature T
// Output: - List with values for epsilon, m, Cv and chi
arma::vec IsingModel::analytical(double T)
{
  double Z = 12 + 2*std::exp(-8/T) + 2*std::exp(8/T); // partition function
  double epsilon = 4/Z * (std::exp(-8/T) - std::exp(8/T)); // epsilon
  double epsilon2 = 8/Z * (std::exp(-8/T) + std::exp(8/T)); // epsilon^2
  double m = std::abs(2/Z * (std::exp(8/T) + 2)); // m
  double m2 = 2/Z * (std::exp(8/T) + 1); // m^2
  double Cv = L_*L_/(T*T) * (epsilon2 - epsilon*epsilon); // heat capacity Cv
  double chi = L_*L_/T * (m2 - m*m); // magn. susceptibility chi
  return arma::vec({epsilon, m, Cv, chi});
}
