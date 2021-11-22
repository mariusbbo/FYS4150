// Markov chain Monte Carlo class

#ifndef __MCMC_hpp__
#define __MCMC_hpp__

#include "armadillo"
#include "IsingModel.hpp"

class MCMC
{
  public:
    int n_Tvals_; // number of temperature values
    arma::vec T_vals_; // list with temperature values
    int n_cycles_; // number of MC cycles
    int L_; // lattice size
    bool ordered_; // determine if initial state is random or ordered
    int n_burnin_; // burn-in time

    // Constructor
    MCMC(double T_start, double T_end, int T_vals, int n_cycles, int L, bool ordered, int n_burnin);

    // Run Markov chain Monte Carlo algorithm
    void run_MCMC();

};

#endif
