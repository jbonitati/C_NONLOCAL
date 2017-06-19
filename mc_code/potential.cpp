#include <cmath>
#include "potential.h"

const double hbarc = 197.3269718;
const double mass_unit=931.494;

/*Returns the potential due to the centrifugal term in the Hamiltonian*/
double Potential::central_potential(double R, Particle proj, Particle targ, Channel * c)
{
  double mu = targ.getM()*proj.getM()/(targ.getM() + proj.getM());
	double constant = (pow(hbarc,2)/(2.0*mass_unit*mu));
  double l = c->getL();

	return constant*((l*(l+1.0))/(R*R));
}

/*Returns the Coulomb potential for particles at distance R*/
double Potential::coulomb_potential(double R, Particle proj, Particle targ)
{
  double rc = 1.25; //don't know how to calculate this
  double Rcoul = rc*pow((targ.getM()),1.0/3.0);
  double z1 = targ.getZ();
  double z2 = proj.getZ();

  if (R<Rcoul)
  {
    return ((z1*z2*1.43997)/(2.0))*(3.0-pow(R,2.0)/(pow(Rcoul,2)));
  }
  else
  {
    return (z1*z2*1.43997)/R;
  }
}  


double V_Potential::getValue(double R, Particle p){
  double Rr = r*pow(p.getM(),1.0/3.0);
  return coeff*-1.0/(1.0+exp((R-Rr)/a));
}

double D_Potential::getValue(double R, Particle p){
  double Rr = r*pow(p.getM(),1.0/3.0);
    return (-4.0*coeff*exp((R-Rr)/a))/pow((1.0+exp((R-Rr)/a)),2.0);
}

double SO_Potential::getValue(double R, Particle p, Channel * c){
  double Rr = r*pow(p.getM(),1.0/3.0);
  double l = c->getL();
  double j = c->getJ();
  double m = c->getSpin();
  
  double Vspin=(j*(j+1.0)-l*(l+1.0)-m*(m+1.0))
  *(-coeff/(a*R))*(pow(hbarc/139.6,2))
  *(exp((R-Rr)/a))/(pow((1.0+exp((R-Rr)/a)),2.0));

  return Vspin;
}
