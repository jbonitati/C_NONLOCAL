#include <cmath>
#include <gsl/gsl_sf.h>
#include <complex>
#include <armadillo>
#include "channel.h"
#include "particle.h"

const double hbarc = 197.3269718; //MeV * fm / c
const double mass_unit=931.454; // MeV / c^2
const double E2HC = 0.007297353; //e^2 / (4pi*epsilon_0) in units of hbar*c
const arma::cx_double I(0.0,1.0);

/*Returns the potential due to the centrifugal term in the Hamiltonian*/
double Channel::central_potential(double R)const
{
	double constant = (pow(hbarc,2)/(2.0*mass_unit*mu));

	return constant*((l*(l+1.0))/(R*R));
}

//returns the wave number of the channel with respect to
// the total energy of the system
//units: fm^{-1}
double Channel::getKc(double energy)const{
  return sqrt(2*mu*mass_unit*std::abs(energy - E)) / hbarc;
}

//returns the relative velocity of the channel with respect to 
// the total energy of the sytem
//units: c (speed of light)
double Channel::getVc(double energy)const{
  double kc = getKc(energy);
  double vc = hbarc*kc / (mu*mass_unit);
  return vc;
}

//returns the Sommerfeld parameter for the channel
//dimensionless
double Channel::getEta(double energy, Particle targ, Particle proj)const{
  double vc = getVc(energy);
  return targ.getZ()*proj.getZ()*E2HC/(hbarc*vc);
}

//calculates the conjugate functions I and O from the coulomb functions at kc*a
//stores the values in I and O, and stores the derivative values in Ip and Op
void Channel::io_coulomb_functions(double x, double energy, Particle targ, Particle proj,
  arma::cx_double *I1, arma::cx_double *O1, arma::cx_double *Ip1, arma::cx_double *Op1)const{

  double etac = getEta(energy,targ,proj);
    
  //double hbarx = hbarc*hbarc/(2*mass_unit*mu);
  //double q = sqrt(energy/hbarx);
  
  gsl_sf_result F1,G1,Fp1,Gp1;
  double exp_F1, exp_G1;
  gsl_sf_coulomb_wave_FG_e(etac,x,l,0,&F1,&Fp1,&G1,&Gp1,&exp_F1,&exp_G1);
  *I1 = G1.val - I*F1.val;
  *Ip1 = (Gp1.val - I*Fp1.val);
  *O1 = G1.val + I*F1.val;
  *Op1 = (Gp1.val + I*Fp1.val);
  //removed factor of "q" from the derivatives
  
}
