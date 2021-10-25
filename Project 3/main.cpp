#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "armadillo"

int main()
{
  // particle parameters
  double x0 = 1;
  double z0 = 1;
  double v0 = 1;
  double q = 1;
  double m = 40;
  arma::vec r = arma::vec({x0, 0, z0});
  arma::vec v = arma::vec({0, v0, 0});

  // Penning trap parameters
  double B0 = 9.65*10; // magnetic field strength
  // double V0 = 9.65*std::pow(10, 8); // old value for potential
  double V0 = 0.0025*9.64852558*1e7; // new value for potential
  // double d = std::pow(10, 4); // old value for char. dimension
  double d = 0.05*1e4; // new value for char. dimension
  double t0 = 0; // intialize time

  double dt = 0.01; // step size
  double t_end = 500; // total simulation time in μm
  int n_timesteps = t_end/dt; // number of timesteps
  int n_particles = 100; // number of particles

  // set equal to true to use RK4 and false to use forward Euler
  bool RK4 = true;

  // set equal to true to turn on Coulomb interactions
  bool Interactions = false;

  // parameters for time-dependent electric field
  // if f=0 the time-independent electric field is used
  double f = 0;
  double omega_v = 0;

  // run simulation with two particles
  // PenningTrap Ptrap = PenningTrap(B0, V0, d, Interactions, f, omega_v, t0);
  // Ptrap.add_particle(Particle(q, m, r, v));
  // Ptrap.add_particle(Particle(q, m, arma::vec({0,1,1}), arma::vec({1,0,0})));
  // Ptrap.simulation(dt, n_particles, n_timesteps, RK4);

  // lists with parameters for time-dependent electric field
  arma::vec f_list = arma::vec({0.1, 0.4, 0.7}); // amplitude list
  arma::vec omega_v_list = arma::regspace(0.4, 0.001, 0.5); // frequency list
  std::cout << "Total number of omega_v's (iterations): " << omega_v_list.n_elem << "\n";

  // vector for storing number of particles inside trap
  arma::mat np_in_trap = arma::mat(3, omega_v_list.n_elem, arma::fill::zeros);

  // run simulation for different values of f
  for (int i = 0; i < f_list.n_elem; i++)
  {
    std::cout << i+1 << "/" << 3 << std::endl;

    // run simulation for different values of omega_v
    for (int j = 0; j < omega_v_list.n_elem; j++)
    {
      // initialize parameters of PenningTrap class
      PenningTrap Ptrap = PenningTrap(B0, V0, d, Interactions, f_list(i), omega_v_list(j), t0);

      // set seed for random number generator
      arma::arma_rng::set_seed(1);

      // generate random initial position and velocity for all particles
      for (int k = 0; k < n_particles; k++)
      {
        r = arma::vec(3).randn() * 0.1 * Ptrap.d_;  // random initial position
        v = arma::vec(3).randn() * 0.1 * Ptrap.d_;  // random initial velocity
        Ptrap.add_particle(Particle(q, m, r, v));
      }
      // run simulation
      Ptrap.simulation(dt, n_particles, n_timesteps, RK4);

      // run function counting number of particles inside trap and store the value in a vector
      np_in_trap(i, j) = Ptrap.count_particles();
      std::cout << "\rIteration " << j+1 << "  " << std::flush;
    }
  }
  // save number of particles inside trap to file
  // np_in_trap.save("particles_in_trap_zoom_ints_test_new.dat");

  // compute analytical solution and save positions and velocities to file
  // Ptrap.analytical_solution(x0, z0, v0, m, q, B0, V0, d, t_end, n_timesteps);

  return 0;
}
