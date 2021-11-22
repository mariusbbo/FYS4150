// Ising model class

#ifndef __IsingModel_hpp__
#define __IsingModel_hpp__

#include "armadillo"

class IsingModel
{
public:
  arma::mat S_; // lattice
  int L_; // lattice size
  arma::vec i_list_; // index_list

  // Constructor
  // Creates initial state
  IsingModel(int L, bool ordered);

  // Compute energy of initial state
  double initialize_energy();

  // Compute magnetization of initial state
  double initialize_magnetization();

  // Create list with indices used for determining probability ratios and change in energy
  void create_index_list();

  // Find the index of the index list used for prob. ratios and change in energy
  int get_index(int i, int j);

  // Create list with probability ratios
  arma::vec probability_ratios(double T);

  // Create list with the possible changes in energy
  arma::vec delta_E();

  // Find position of the current energy value when creating a histogram
  int find_histogram_pos(arma::vec comp_list, arma::vec test, double E);

  // Compute analytical values
  arma::vec analytical(double T);

};

#endif
