#include <complex>
#include "potential.h"

class OpticalPotential{
  private:
    const V_Potential Volume_Re, Volume_Im;
    const D_Potential Surface_Re, Surface_Im;
    const SO_Potential SpinOrbit_Re, SpinOrbit_Im;
      
  protected:
  
    std::complex<double> Volume_potential(double radius, Particle p);
    std::complex<double> Surface_potential(double radius, Particle p);
    std::complex<double> Spin_orbit_potential(double radius, Particle p, Channel * c);
    
  public:
    OpticalPotential(const V_Potential Vr, const V_Potential Vi, 
      const D_Potential Sr, const D_Potential Si, 
      const SO_Potential SOr, const SO_Potential SOi):
      Volume_Re(Vr), Volume_Im(Vi), Surface_Re(Sr), Surface_Re(Si),
      SpinOrbit_Re(SOr), SpinOrbit_Im(SOi){ }
      
    //todo: write explicit copy constructor
      
    std::complex<double> totalPotential(double radius, Particle p, Channel * c);
};

class NonLocalOpticalPotential : public OpticalPotential{
  private:
    const double beta;//range of nonlocality

  public:
    NonLocalOpticalPotential(const Potential Vr, const Potential Vi, 
      const Potential Sr, const Potential Si, 
      const Potential SOr, const Potential SOi, const double beta_):
      OpticalPotential(Vr,Vi,Sr,Si,SOr,SOi), beta(beta_){ }
      
    std::complex<double> totalPotential(double r1, double r2, Particle p, Channel * c);
};
