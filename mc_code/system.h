#include "channel.h"
#include "particle.h"
#include "OpticalPotential.h"
#include "lagbasis.h"
#include "wavefunction.h"
#include "potential.h"
#include "constants.h"

#include <armadillo>
#include <vector>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/path.hpp>
//#include <boost/serialization/array_wrapper.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas; //matrix

class System{
  private:
    const double a; // channel radius
    double energy; //projectile energy
    
    const Particle proj;  //projectile particle
    const Particle targ; //target particle
    
    double coulomb_radius; //Coulomb radius (rc)
    
    const unsigned int entrance_channel;
    std::vector<Channel> channels; //an array of channels
    
    const OpticalPotential pot; //manages local potential
    const NonLocalOpticalPotential nlpot; //manages nonlocal potential
    
    int basis_size;
    double step_size;
    double r_max;
    LagBasis lagrangeBasis;
    
    double beta; //coupling constant
    
    arma::cx_mat cmatrix; //C and its inverse are matrices of matrices
    arma::cx_mat invcmatrix;
    arma::cx_mat rmatrix;
    arma::cx_mat umatrix;
    
    
    matrix<CouplingPotential> coupling_matrix;
    
    void cmatrixCalc();
    void invertCmatrix();
    void rmatrixCalc();
    void umatrixCalc();
    
    double couplingPotential(double, double, unsigned int, unsigned int);
    
    WaveFunction calculateWaveFunction(double);
    
    System(); //explicitly disallow unwanted default constructor (Meyers item 6)
  public:
    ~System(){
      channels.clear();
    }
    
    arma::cx_mat getCmatrix(){return cmatrix;}
    arma::cx_mat getInvCmatrix(){return invcmatrix;}
    arma::cx_mat getRmatrix(){return rmatrix;}
    arma::cx_mat getUmatrix(){return umatrix;}
    
    System(const double a_size, double e,
      Particle projectile, Particle target, double rc,
      int c0, std::vector<Channel> channels_,
      OpticalPotential op, NonLocalOpticalPotential nlop, int basis_size,
      double step, double max,matrix<CouplingPotential> coupling);
    
    cx_double internalWaveFunction(unsigned int, double);
    cx_double externalWaveFunction(unsigned int, double);
    
    void waveFunctions(boost::filesystem::ofstream&);
    void plotPotential(boost::filesystem::ofstream&);
};
