#include <armadillo>
#include "potential.h"
#include "particle.h"

class OpticalPotential{
  private:
    const V_Potential Volume_Re, Volume_Im;
    const D_Potential Surface_Re, Surface_Im;
    const SO_Potential SpinOrbit_Re, SpinOrbit_Im;
      
  protected:
  
    arma::cx_double Volume_potential(double radius, Particle p)const;
    arma::cx_double Surface_potential(double radius, Particle p)const;
    arma::cx_double Spin_orbit_potential(double radius, Particle p, Channel * c)const;
    
  public:
    OpticalPotential(const Potential Vr, const Potential Vi, 
      const Potential Sr, const Potential Si, 
      const Potential SOr, const Potential SOi):
      Volume_Re(Vr), Volume_Im(Vi), Surface_Re(Sr), Surface_Im(Si),
      SpinOrbit_Re(SOr), SpinOrbit_Im(SOi){ }
      
    //todo: write explicit copy constructor
      
    arma::cx_double totalPotential(double radius, Particle p, Channel * c)const;
};

class NonLocalOpticalPotential : public OpticalPotential{
  private:
    const double beta;//range of nonlocality

  public:
    NonLocalOpticalPotential(const Potential Vr, const Potential Vi, 
      const Potential Sr, const Potential Si, 
      const Potential SOr, const Potential SOi, const double beta_):
      OpticalPotential(Vr,Vi,Sr,Si,SOr,SOi), beta(beta_){ }
      
    arma::cx_double totalPotential(double r1, double r2, Particle p, Channel * c)const;
};
