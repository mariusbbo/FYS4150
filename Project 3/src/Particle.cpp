#include "Particle.hpp"

// Constructor
Particle::Particle(double charge, double mass, arma::vec position, arma::vec velocity)
{
  q_ = charge;
  m_ = mass;
  r_ = position;
  v_ = velocity;
}
