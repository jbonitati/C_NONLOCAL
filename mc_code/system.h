#include "channel.h"
#include "particle.h"
#include "OpticalPotential.h"
#include "lagbasis.h"

#include <vector>
#include <armadillo>
#include <complex>
#include <cstring>


class System{
  private:
    const double a; // channel radius
    double energy; //projectile energy
    
    const Particle proj;  //projectile particle
    const Particle targ; //target particle
    
    const int num_channels;
    std::vector<Channel*> channels;
    
    const OpticalPotential pot; //manages local potential
    const NonLocalOpticalPotential nlpot; //manages nonlocal potential
    
    int basis_size;
    double step_size;
    double r_max;
    LagBasis lagrangeBasis;
    
    double beta; //coupling constant
    
    arma::mat<arma::cx_mat> cmatrix;
    arma::mat<arma::cx_mat> invcmatrix;
    arma::cx_mat rmatrix;
    arma::cx_mat umatrix;
    
    void cmatrixCalc();
    void rmatrixCalc();
    void umatrixCalc();
    
    System(){ } //explicitly disallow default constructor
  public:
    ~System(){
      for (std::vector<Channel *>::iterator it = channels.begin() ; it != channels.end(); ++it)
        delete (*it);
        
      channels.clear();
    }
    
    System(const double a_size, double e,
      const double m1, const double m2, const double z1, const double z2,
      const int num_channels_, std::vector<Channel*> channels_,
      OpticalPotential op, NonLocalOpticalPotential nlop, int basis_size,
      double step, double max, double coupling);
    
    void waveFunction(ofstream outFile);
};
