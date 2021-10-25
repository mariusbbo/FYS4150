#include "PenningTrap.hpp"
#include "Particle.hpp"
#include "math.h"

// Constructor
PenningTrap::PenningTrap(double B0, double V0, double d, bool Interactions, double f, double omega_v, double t)
{
  B0_ = B0;
  V0_ = V0;
  d_ = d;
  ParticleObjects_ = std::vector<Particle> {};
  Interactions_ = Interactions;
  f_ = f;
  omega_v_ = omega_v;
  t_ = t;
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in)
{
  // add particle objects, charge, mass, position and velocity, to list
  ParticleObjects_.push_back(p_in);
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r)
{
  // check if particle is outside trap
  if (arma::norm(r) > d_)
  {
    // set electric field to zero outside trap
    return arma::vec({0,0,0});
  }
  else
  {
  double V = V0_ * (1 + f_*cos(omega_v_*t_)); // potential
  double C = V / std::pow(d_,2); // V0/d^2
  double Ex = C * r(0);
  double Ey = C * r(1);
  double Ez = -C * 2 * r(2);
  arma::vec E = {Ex, Ey, Ez}; // electric field vector
  return E;
  }
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r)
{
  // check if particle is outside trap
  if (arma::norm(r) > d_)
  {
    // set magnetic field to zero outside trap
    return arma::vec({0,0,0});
  }
  else
  {
  arma::vec B = {0, 0, B0_}; // magnetic field vector
  return B;
  }
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j)
{
  double ke = 1.38935333 * std::pow(10,5); // Coulomb constant
  double qi = ParticleObjects_[i].q_; // charge of particle_i
  double qj = ParticleObjects_[j].q_; // charge of particle_j
  double mi = ParticleObjects_[i].m_; // mass of particle_i

  arma::vec ri = ParticleObjects_[i].r_; // position of particle_i
  arma::vec rj = ParticleObjects_[j].r_; // position of particle_j
  arma::vec r_diff = ri - rj;

  // find distance between particle_i and particle_j
  double rx = std::pow(r_diff(0), 2);
  double ry = std::pow(r_diff(1), 2);
  double rz = std::pow(r_diff(2), 2);
  double r = std::pow((rx+ry+rz), 0.5);

  // compute components of Coulomb force between particle_i and particle_j
  double fx = ke * qi*qj / mi * (ri(0) - rj(0)) / std::pow(r, 3);
  double fy = ke * qi*qj / mi * (ri(1) - rj(1)) / std::pow(r, 3);
  double fz = ke * qi*qj / mi * (ri(2) - rj(2)) / std::pow(r, 3);
  arma::vec Fp = {fx, fy, fz}; // Coulomb force vector

  return Fp;
}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i)
{
  double q = ParticleObjects_[i].q_; // charge
  double m = ParticleObjects_[i].m_; // mass
  arma::vec r = ParticleObjects_[i].r_; // position
  arma::vec v = ParticleObjects_[i].v_; // velocity

  arma::vec E = external_E_field(r); // compute electric field at position r
  arma::vec B = external_B_field(r); // compute magnetic field at position r

  // compute Lorentz force
  arma::vec F_L = (q*E + q*arma::cross(v, B));
  return F_L;
}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i)
{
  int N = ParticleObjects_.size(); // number of particles
  arma::vec Fp_tot = arma::vec(3, arma::fill::zeros);
  for (int j = 0; j < N; j++)
  {
    if (i != j) // only compute force on particle_i from the other particles
    {
      // sum up all forces from the other particles
      Fp_tot += force_particle(i, j);
    }
  }
  return Fp_tot;

}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i)
{
  // check which determines if simulation is run with or without particle interactions
  if (Interactions_ == true)
  {
    return total_force_external(i) + total_force_particles(i);
  }
  else
  {
    return total_force_external(i);
  }
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt)
{
  double t = t_; // temporary time that will not be changed during steps
  double m; // mass
  arma::vec r; // position
  arma::vec v; // velocity
  arma::vec k1r, k1v, k2r, k2v, k3r, k3v, k4r, k4v; // intermediate positions and velocities for RK4
  int N = ParticleObjects_.size(); // number of particles

  // evolve all particles one time step
  for (int i = 0; i < N; i++)
  {
    m = ParticleObjects_[i].m_; // mass
    r = ParticleObjects_[i].r_; // position
    v = ParticleObjects_[i].v_; // velocity

    // compute k1, k2, k3 and k4 for both velocity and position and
    // update position and velocity and time for each intermediate step

    k1v = dt * 1/m * total_force(i);
    k1r = dt * v;

    ParticleObjects_[i].v_ = v + 0.5 * k1v;
    ParticleObjects_[i].r_ = r + 0.5 * k1r;
    t_ = t + 0.5*dt;

    k2v = dt * 1/m * total_force(i);
    k2r = dt * (v + 0.5*k1v);

    ParticleObjects_[i].v_ = v + 0.5 * k2v;
    ParticleObjects_[i].r_ = r + 0.5 * k2r;
    t_ = t + 0.5*dt;

    k3v = dt * 1/m * total_force(i);
    k3r = dt * (v + 0.5*k2v);

    ParticleObjects_[i].v_ = v + k3v;
    ParticleObjects_[i].r_ = r + k3r;
    t_ = t + dt;

    k4v = dt * 1/m * total_force(i);
    k4r = dt * (v + k3v);

    ParticleObjects_[i].v_ = v + 1/6.*(k1v + 2*k2v + 2*k3v + k4v); // update position for time step i+1
    ParticleObjects_[i].r_ = r + 1/6.*(k1r + 2*k2r + 2*k3r + k4r); // update velocity for time step i+1
  }
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{
  double m; // mass
  int N = ParticleObjects_.size(); // number of particles

  // evolve all particles one time step
  for (int i = 0; i < N; i++)
  {
    m = ParticleObjects_[i].m_;
    ParticleObjects_[i].v_ = ParticleObjects_[i].v_ + dt* 1/m*total_force(i);
    ParticleObjects_[i].r_ = ParticleObjects_[i].r_ + dt*ParticleObjects_[i].v_;
  }
}

// Simulates Penning trap with Runge Kutta 4 or Forward Euler and save positions and velocities to file
void PenningTrap::simulation(double dt, int n_particles, int n_timesteps, bool RK4)
{
  // vectors for storing positions and velocities
  arma::cube position = arma::cube(3, n_particles, n_timesteps, arma::fill::zeros);
  arma::cube velocity = arma::cube(3, n_particles, n_timesteps, arma::fill::zeros);

  // run Runge Kutta 4 or Forward Euler for all time steps
  for (int i = 0; i < n_timesteps; i++)
  {
    // check if Runge Kutta 4 or Forward Euler is be used
    if (RK4 == true)
    {
      evolve_RK4(dt);
    }
    else
    {
      evolve_forward_Euler(dt);
    }
    for (int j = 0; j < n_particles; j++)
    {
      // add position and velocity for all particles at each time step
      position.slice(i).col(j) = ParticleObjects_[j].r_;
      velocity.slice(i).col(j) = ParticleObjects_[j].v_;
    }
  }
  // save positions and velocities to files
  // position.save("position.dat");
  // velocity.save("velocity.dat");
}

// Compute analytical solution and save positions to file
void PenningTrap::analytical_solution(double x0, double z0, double v0, double m,
                                      double q, double B0, double V0, double d,
                                      double t_end, int n_timesteps)
{
  double omega0 = q*B0/m;
  double omega_z = 2*q*V0/(m*std::pow(d,2));

  double omega_p = (omega0 + std::sqrt(std::pow(omega0,2) - 2*omega_z)) / 2; // omega_+
  double omega_m = (omega0 - std::sqrt(std::pow(omega0,2) - 2*omega_z)) / 2; // omega_-

  double Ap = (v0 + omega_m*x0) / (omega_m - omega_p); // A_+
  double Am = -(v0 + omega_p*x0) / (omega_m - omega_p); // A_-

  arma::vec t = arma::linspace(0, t_end, n_timesteps);

  arma::vec x = Ap*cos(omega_p*t) + Am*cos(omega_m*t);
  arma::vec y = -Ap*sin(omega_p*t) - Am*sin(omega_m*t);
  arma::vec z = z0*cos(std::sqrt(omega_z)*t);

  arma::mat r = arma::mat(3, n_timesteps, arma::fill::zeros);

  // add solutions for x(t), y(t) and z(t) to vector
  for (int i = 0; i < n_timesteps; i++)
  {
    r.col(i) = arma::vec({x(i), y(i), z(i)});
  }
  // save analytical solution to file
  r.save("position_analytical.dat");
}

// Count the number of particles remaining inside the trap
int PenningTrap::count_particles()
{
  int n_in = 0; // counter
  int N = ParticleObjects_.size(); // total number of particles

  // loop through positions of all particles
  for (int i = 0; i < N; i++)
  {
    // add 1 to counter if the particle is inside the trap
    if (arma::norm(ParticleObjects_[i].r_) < d_)
    {
      n_in += 1;
    }
  }
  return n_in;
}
