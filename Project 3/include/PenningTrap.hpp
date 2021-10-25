// Penning trap class

#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include "armadillo"
#include "Particle.hpp"

class PenningTrap
{
  public:
    double B0_; // magnetic field strength
    double V0_; // applied potential
    double d_; // characteristic dimension
    std::vector<Particle> ParticleObjects_; // vector for storing all particle objects
    bool Interactions_; // decides if Coulomb interactions between particles occur
    double f_; // amplitude of time-dependent potential
    double omega_v_; // frequency for time-dependent potential
    double t_; // initialize time

    // Constructor
    PenningTrap(double B0, double V0, double d, bool Interactions, double f, double omega_v, double t);

    // Add a particle to the trap
    void add_particle(Particle p_in);

    // External electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r);

    // External magnetic field at point r=(x,y,z)
    arma::vec external_B_field(arma::vec r);

    // Force on particle_i from particle_j
    arma::vec force_particle(int i, int j);

    // The total force on particle_i from the external fields
    arma::vec total_force_external(int i);

    // The total force on particle_i from the other particles
    arma::vec total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    arma::vec total_force(int i);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt);

    // Run evolve_forward_Euler or evolve_RK4
    void simulation(double dt, int n_particles, int n_timesteps, bool RK4);

    // Compute analytical solution
    void analytical_solution(double x0, double z0, double v0, double m,
                                          double q, double B0, double V0, double d,
                                          double t_end, int n_timesteps);

    // Count the number of particles inside the trap
    int count_particles();

};

#endif
