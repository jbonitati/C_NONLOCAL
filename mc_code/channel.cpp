#include <cmath>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_result.h>
#include <armadillo>
#include "channel.h"
#include "particle.h"
#include "constants.h"

//Channel Constructor
//uses spin = total_angular_momentum - angular_momentum
//sets the bloch constant to 0 for open channels or the logarithmic derivative
// of the whittaker function for closed channels
Channel::Channel(const double ec, const int l_, const double m_, const double j_,
  const double mu_, double energy_, double a, Particle targ, Particle proj):
    E(ec), l(l_), j(j_), m(m_), mu(mu_), b(0), energy(energy_){ 
    if(energy < ec){
      double kc = getKc();
      b = 2*kc*a*whittaker_prime(2*kc*a, targ, proj)/whittaker(2*kc*a, targ, proj);
    }
}

/*Returns the potential due to the centrifugal term in the Hamiltonian*/
double Channel::central_potential(double R)const{
	return (pow(hbarc,2)/(2.0*mass_unit*mu))*((l*(l+1.0))/(R*R));
}

//returns the wave number of the channel with respect to
// the total energy of the system
//units: fm^{-1}
double Channel::getKc()const{
  return sqrt(2*mu*mass_unit*std::abs(energy-E))
    / hbarc;
}

//returns the relative velocity of the channel with respect to 
// the total energy of the sytem
//units: c (speed of light)
double Channel::getVc()const{
  double kc = getKc();
  double vc = hbarc*kc / (mu*mass_unit);
  return vc;
}

//returns the Sommerfeld parameter for the channel
//dimensionless
double Channel::getEta(Particle targ, Particle proj)const{
  double vc = getVc();
  return targ.getZ()*proj.getZ()*E2HC/(hbarc*vc);
}

//calculates the conjugate functions I and O from the coulomb functions at kc*a
//stores the values in I and O, and stores the derivative values in Ip and Op
void Channel::io_coulomb_functions(double x, Particle targ, Particle proj,
  arma::cx_double *I1, arma::cx_double *O1, arma::cx_double *Ip1, arma::cx_double *Op1)const{

  double etac = getEta(targ,proj);
  
  gsl_sf_result F1,G1,Fp1,Gp1;
  double exp_F1, exp_G1;
  gsl_sf_coulomb_wave_FG_e(etac,x,l,0,&F1,&Fp1,&G1,&Gp1,&exp_F1,&exp_G1);
  double G = G1.val*exp(exp_G1);
  double F = F1.val*exp(exp_F1);
  double Gp = Gp1.val*exp(exp_G1);
  double Fp = Fp1.val*exp(exp_F1);
  *I1 = G - I*F;
  *Ip1 = Gp - I*Fp;
  *O1 = G + I*F;
  *Op1 = Gp + I*Fp;
  //removed a factor of k from the derivatives that shouldn't have been there
  
}

/*function to test the gsl coulomb functions*/
void Channel::coulomb_test(double x, Particle targ, Particle proj)const{
  double etac = getEta(targ,proj);
  gsl_sf_result F1,G1,Fp1,Gp1;
  double exp_F1, exp_G1;
  gsl_sf_coulomb_wave_FG_e(etac,x,l,0,&F1,&Fp1,&G1,&Gp1,&exp_F1,&exp_G1);
  double G = G1.val*exp(exp_G1);
  double F = F1.val*exp(exp_F1);
  double Gp = Gp1.val*exp(exp_G1);
  double Fp = Fp1.val*exp(exp_F1);
  std::cout << "G " << G << "\nF " << F <<  "\nG' " << Gp << "\nF' " << Fp << std::endl;
  std::cout << exp_F1 << " " << exp_G1 << std::endl;
}

double Channel::whittaker(double k, double m, double z){
  double coeff = exp(-1.0*z/2)*pow(z,m+0.5);
  return coeff * gsl_sf_hyperg_U(0.5+m-k, 1+2*m, z);
}

//whittaker function
//used in calculating the wave function and bloch operator
double Channel::whittaker(double z, Particle targ, Particle proj){
  double k, m;
  k = -1.0*getEta(targ, proj);
  m = l + 0.5;
  return whittaker(k,m,z);
}

//derivative of the whittaker function
//used in calculating the bloch operator for closed channels
double Channel::whittaker_prime(double z, Particle targ, Particle proj){
  double k, m;
  k = -1.0*getEta(targ, proj);
  m = l + 0.5;
  return ((z-2*k)*whittaker(k,m,z) - 2*whittaker(k+1,m,z))/(2*z);
}

