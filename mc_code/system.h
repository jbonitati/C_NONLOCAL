#include "channel.h"
#include "particle.h"
#include "OpticalPotential.h"
#include "lagbasis.h"

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
    
    matrix<Coupling_Potential> coupling_matrix;
    
    void cmatrixCalc();
    void rmatrixCalc();
    void umatrixCalc();
    
    double couplingPotential(double r1, double r2, int c1, int c2);
    
    System(); //explicitly disallow default constructor
  public:
    ~System(){
      channels.clear();
    }
    
    System(const double a_size, double e,
      const double m1, const double m2, const double z1, const double z2,
      int c0, std::vector<Channel> channels_,
      OpticalPotential op, NonLocalOpticalPotential nlop, int basis_size,
      double step, double max, matrix<Coupling_Potential> coupling);
    
    void waveFunction(boost::filesystem::ofstream& outFile);
};
