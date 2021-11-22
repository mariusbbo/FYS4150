#include "IsingModel.hpp"
#include "MCMC.hpp"
#include "armadillo"
#include "chrono"

int main()
{

  // arma::vec n_list({100, 1000, 1e4, 1e5, 1e6}); // list with values for number of MC cycles
  // arma::mat params(2, 5, arma::fill::zeros); // list for storing values for timing tests
  //
  // // loop for running through values in n_list for timing tests
  // for (int i = 0; i < 5; i++)
  // {
  //   // Start measuring time
  //   auto t1 = std::chrono::high_resolution_clock::now();

    double T_start = 2.1; // lowest temperature
    double T_end = 2.4; // highest temperature
    double n_Tvals = 2; // number of temperature values
    int n_cycles = 5e5; // number of MC cycles
    int L = 80; // lattice size
    bool ordered = false; // set to true to initialize ordered state
    int n_burnin = 0; // number of cycles that will be skipped due to burn-in

    // set parameters for Markov chain Monte Carlo class
    MCMC Mc = MCMC(T_start, T_end, n_Tvals, n_cycles, L, ordered, n_burnin);

    // run Markov chain Monte Carlo algorithm n_cycles times for n_Tvals tempeature values
    Mc.run_MCMC();

  //   // Stop measuring time
  //   auto t2 = std::chrono::high_resolution_clock::now();
  //
  //   // Calculate the elapsed time
  //   // We use chrono::duration<double>::count(), which by default returns duration in seconds
  //   double duration_seconds = std::chrono::duration<double>(t2 - t1).count();
  //   std::cout << "time " << duration_seconds << "\n";
  //   params(0,i) = n_list(i);
  //   params(1,i) = duration_seconds;
  // }
  // params.save("n_time_list_para.dat");


  return 0;
}
