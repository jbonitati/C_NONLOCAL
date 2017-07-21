#include "channel.h"
#include "particle.h"
#include "OpticalPotential.h"
#include "lagbasis.h"
#include "wavefunction.h"
#include "potential.h"
#include "constants.h"

#include <boost/serialization/array_wrapper.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <fstream>
#include <armadillo>
#include <complex>
#include <cstring>
#include <vector>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/path.hpp>

using namespace boost::numeric::ublas;

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
    WaveFunction internalWaveFunction(double);
    WaveFunction externalWaveFunction(double);
    
    System(); //explicitly disallow unwanted default constructor (Meyers item 6)
  public:
    ~System(){
      channels.clear();
    }
    
    System(const double a_size, double e,
      Particle projectile, Particle target, double rc,
      int c0, std::vector<Channel> channels_,
      OpticalPotential op, NonLocalOpticalPotential nlop, int basis_size,
      double step, double max,matrix<CouplingPotential> coupling);
    
    void waveFunctions(boost::filesystem::ofstream&);
    void plotPotential(boost::filesystem::ofstream&);
};
