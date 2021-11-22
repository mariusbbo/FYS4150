#include "IsingModel.hpp"
#include "MCMC.hpp"
#include "armadillo"
#include "random"
#include "omp.h"

// constructor
MCMC::MCMC(double T_start, double T_end, int n_Tvals, int n_cycles, int L, bool ordered, int n_burnin)
{
  n_Tvals_ = n_Tvals; // number of temperature values
  n_cycles_ = n_cycles; // number of MC cycles
  T_vals_ = arma::linspace(T_start, T_end, n_Tvals_); // temperature list
  L_ = L; // lattice size
  ordered_ = ordered; // set initial state to be random or ordered
  n_burnin_ = n_burnin; // burn-in time
}

// Run the Markov chain Monte Carlo algorithm
// Sample states and uses the IsingModel class to compute energy,
// magnetization, heat capacity and magnetic susceptibilty
void MCMC::run_MCMC()
{
  double N = L_*L_; // number of spins in lattice
  // lists used to study burn-in time
  arma::mat epsilon(n_Tvals_, n_cycles_, arma::fill::zeros); // epsilon list
  arma::mat m(n_Tvals_, n_cycles_, arma::fill::zeros); // m list

  // lists storing analytical and numerical results. Store <epsilon>, <m>, Cv and chi for all temperatures
  arma::mat analytical_results(n_Tvals_, 4, arma::fill::zeros);
  arma::mat numerical_results(n_Tvals_, 4, arma::fill::zeros);

  // values for creating histogram
  int N_bins = 800; // number of bins
  double max_eps = 4*N; // maximum energy
  double min_eps = -4*N; // minimum energy
  arma::mat E_histogram(n_Tvals_, N_bins, arma::fill::zeros); // histogram list for all temperatures
  arma::vec E_interval = arma::linspace(min_eps, max_eps, N_bins); // energy interval for histogram
  arma::vec test(N_bins, arma::fill::ones); // list with ones

  // start parallelization
  #pragma omp parallel
  {
  double A; // value for acceptance rule
  double E; // energy
  double M; // magnetization
  double E_temp; // temporary energy
  double M_temp; // temporary magnetization
  arma::vec E_list(n_cycles_, arma::fill::zeros); // energy list
  arma::vec E2_list(n_cycles_, arma::fill::zeros); // list for energy squared
  arma::vec M_list(n_cycles_, arma::fill::zeros); // magnetization list
  arma::vec M2_list(n_cycles_, arma::fill::zeros); // list for magnetization squared


  // define random number generator
  int thread_number = omp_get_thread_num();
  std::mt19937 gen(17199+thread_number);
  std::uniform_real_distribution<double> distribution(0, 1);

  #pragma omp for
  // loop over all temperature values
  for (int i = 0; i < n_Tvals_; i++)
  {
    std::cout << i << std::endl;
    // create instance of IsingModel class
    IsingModel I_model = IsingModel(L_, ordered_);
    arma::vec deltaE = I_model.delta_E(); // create list with possible changes in energy

    // compute analytical results
    analytical_results(i,0) = I_model.analytical(T_vals_(i))(0);
    analytical_results(i,1) = I_model.analytical(T_vals_(i))(1);
    analytical_results(i,2) = I_model.analytical(T_vals_(i))(2);
    analytical_results(i,3) = I_model.analytical(T_vals_(i))(3);

    I_model.create_index_list(); // create index list for prob. ratios and change in energy
    arma::vec p_ratio_list = I_model.probability_ratios(T_vals_(i)); // create list with prob. ratios

    // find energy and magnetization for initial state
    E = I_model.initialize_energy();
    M = I_model.initialize_magnetization();

    // loop over all MC cycles
    for (int j = 0; j < n_cycles_; j++)
    {
      // loop over all spin flips
      for (int k = 0; k < N; k++)
      {
        // generate random position in lattice
        int a = floor(distribution(gen)*L_);
        int b = floor(distribution(gen)*L_);

        int index = I_model.get_index(a,b); // find index for the spin at position (a,b)
        double p_ratio = p_ratio_list(index); // find prob. ratio
        A = std::min(p_ratio, 1.0); // check if prob. ratio is below one
        double r = distribution(gen); // generate random number from uniform dist between 0 and 1

        // acceptance rule, flip spin if A larger than r
        if (A >= r)
        {
          I_model.S_(a,b) *= -1; // flip spin
          E += deltaE(index); // update energy
          M += 2*I_model.S_(a,b); // update magnetization
        }

      }
      // find epsilon and m for study of burn-in time
      E_temp += E;
      M_temp += std::abs(M);
      epsilon(i,j) = E_temp/N/(j+1);
      m(i,j) = M_temp/N/(j+1);

      E_list(j) = E; // energy list
      M_list(j) = std::abs(M); // magnetization list
      E2_list(j) = E*E; // list with energy squared
      M2_list(j) = M*M; // list with magnetization squared

      // create histogram
      int hist_index = I_model.find_histogram_pos(E_interval, test, E); // find position of current energy in histogram
      E_histogram(i, hist_index) += 1; // add 1 to position of current energy

    }
    // sum over E, M, E^2, M^2 for all MC cycles for each temperature
    double E_sum = arma::sum(E_list.rows(n_burnin_, n_cycles_-1));
    double M_sum = arma::sum(M_list.rows(n_burnin_, n_cycles_-1));
    double E2_sum = arma::sum(E2_list.rows(n_burnin_, n_cycles_-1));
    double M2_sum = arma::sum(M2_list.rows(n_burnin_, n_cycles_-1));

    // compute <epsilon>, <|m|>, Cv and chi for each temperature and store in array
    numerical_results(i,0) = E_sum/N/(n_cycles_-n_burnin_); // average energy per spin <epsilon>
    numerical_results(i,1) = M_sum/N/(n_cycles_-n_burnin_); // average magnetization per spin <|m|>
    numerical_results(i,2) = 1/N * 1/(T_vals_(i)*T_vals_(i)) * (E2_sum/(n_cycles_-n_burnin_) - ( std::pow(E_sum/(n_cycles_-n_burnin_), 2)) ); // heat capacity Cv
    numerical_results(i,3) = 1/(N*T_vals_(i)) * (M2_sum/(n_cycles_-n_burnin_) - ( std::pow(M_sum/(n_cycles_-n_burnin_), 2) )); // magnetic susceptibility chi

  }
}

  // save values of <epsilon>, <|m|>, Cv and chi for all temperatures for different lattices sizes
  // numerical_results.save("params_L40.dat");
  // numerical_results.save("params_L60.dat");
  // numerical_results.save("params_L80.dat");
  // numerical_results.save("params_L100.dat");

  // compute and save analytical and numerical results for 2x2 matrix
  // std::cout << "Analytical " << "\n";
  // std::cout << analytical_results << "\n";
  // analytical_results.save("analytical_results.dat");
  //
  // std::cout << "Numerical " << "\n";
  // std::cout << numerical_results << "\n";
  // numerical_results.save("numerical_n1e4.dat");

  // save values for estimating burn-in time
  // epsilon.save("epsilon.dat");
  // m.save("m.dat");
  epsilon.save("epsilon_ordered.dat");
  m.save("m_ordered.dat");

  // save values for histogram
  // E_histogram.save("E_histogram.dat");

}
