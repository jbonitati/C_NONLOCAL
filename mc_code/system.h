#include "channel.h"
#include "particle.h"
#include "OpticalPotential.h"
#include "lagbasis.h"

#include <boost/serialization/array_wrapper.hpp>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <armadillo>
#include <complex>
#include <cstring>
#include <fstream>

using namespace boost::numeric::ublas;

class System{
  private:
    const double a; // channel radius
    double energy; //projectile energy
    
    const Particle proj;  //projectile particle
    const Particle targ; //target particle
    
    const int num_channels;
    int entrance_channel;
    std::vector<Channel*> channels;
    
    const OpticalPotential pot; //manages local potential
    const NonLocalOpticalPotential nlpot; //manages nonlocal potential
    
    int basis_size;
    double step_size;
    double r_max;
    LagBasis lagrangeBasis;
    
    double beta; //coupling constant
    
    matrix<arma::cx_mat> cmatrix; //C and its inverse are matrices of matrices
    matrix<arma::cx_mat> invcmatrix;
    arma::cx_mat rmatrix;
    arma::cx_mat umatrix;
    
    void cmatrixCalc();
    void rmatrixCalc();
    void umatrixCalc();
    
    arma::cx_double couplingPotential(double r1, double r2, Channel * c1, Channel * c2);
    
    System(); //explicitly disallow default constructor
  public:
    ~System(){
      for (std::vector<Channel *>::iterator it = channels.begin() ; it != channels.end(); ++it)
        delete (*it);
      channels.clear();
    }
    
    System(const double a_size, double e,
      const double m1, const double m2, const double z1, const double z2,
      const int num_channels_, int c0, std::vector<Channel*> channels_,
      OpticalPotential op, NonLocalOpticalPotential nlop, int basis_size,
      double step, double max, double coupling);
    
    void waveFunction(std::ofstream& outFile);
};
