#include <cmath>
#include <gsl/gsl_sf>
#include <complex>
#include <armadillo>
#include "channel.h"
#include "particle.h"

const double hbarc = 197.3269718;
const double mass_unit=931.494;
const double E2HC = 0.00729927;
const arma::cx_double I(0.0,1.0);

/*Returns the potential due to the centrifugal term in the Hamiltonian*/
double Channel::central_potential(double R)
{
	double constant = (pow(hbarc,2)/(2.0*mass_unit*mu));

	return constant*((l*(l+1.0))/(R*R));
}

//returns the wave number of the channel with respect to
// the total energy of the system
double Channel::getKc(double energy){
  if(energy > E)
   return sqrt(2*mu*mass_unit*(energy - E)) / hbarc;
  else
   return sqrt(2*mu*mass_unit*(E - energy)) / hbarc;
}

//returns the relative velocity of the channel with respect to 
// the total energy of the sytem
double Channel::getVc(double energy){
  double kc = getKc(energy);
  double vc = hbarc*kc / (mu*mass_unit);
  return vc;
}

//calculates the conjugate functions I and O from the coulomb functions at kc*a
//stores the values in I and O, and stores the derivative values in Ip and Op
void Channel::io_coulomb_functions(double x, double energy, Patricle targ, Particle Proj,
  arma::cx_double *I1, arma::cx_double *O1, arma::cx_double *Ip1, arma::cx_double *Op1){

  double kc = getKc(energy);
  double vc = getVc(energy);
  double etac = targ.getZ()*proj.getZ()*E2HC/(hbarc*vc);
    
  double hbarx = hbarc*hbarc/(2*mass_unit*mu);
  double q = sqrt(energy/hbarx);
  
  gsl_sf_result F1,G1,Fp1,Gp1;
  double exp_F1, exp_G1, exp_F2, exp_G2;
  gsl_sf_coulomb_wave_FG_e(etac,x,L,0,&F1,&Fp1,&G1,&Gp1,&exp_F1,&exp_G1);
  I1 = G1.val - I*F1.val;
  Ip1 = q*(Gp1.val - I*Fp1.val);
  O1 = G1.val + I*Fl.val;
  Op1 = q*(Gp1.val + I*Fp1.val);
  
}

