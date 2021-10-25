// Particle class

#ifndef __Particle_hpp__
#define __Particle_hpp__

#include "armadillo"

class Particle
{
  public:
    double q_; // charge
    double m_; // mass
    arma::vec r_; // position
    arma::vec v_; // velocity

    // Constructor
    Particle(double charge, double mass, arma::vec position, arma::vec velocity);
};

#endif
