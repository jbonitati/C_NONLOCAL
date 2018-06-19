#include <cmath>
#include "potential.h"
#include "constants.h"

double V_Potential::getValue(double R, Particle p)const
{
  double Rr = r;//*pow(p.getM(),1.0/3.0);
  return (coeff == 0 ? 0 : coeff*-1.0/(1.0+exp((R-Rr)/a)));
}

double D_Potential::getValue(double R, Particle p)const
{
  double Rr = r;//*pow(p.getM(),1.0/3.0);
  return (coeff == 0 ? 0 : (-4.0*coeff*exp((R-Rr)/a))/pow((1.0+exp((R-Rr)/a)),2.0));
}

double SO_Potential::getValue(double R, Particle p, Channel * c)const
{
  double Rr = r;//*pow(p.getM(),1.0/3.0);
  int l = c->getL();
  double j = c->getJ();
  double m = c->getSpin();
  
  double Vspin= (coeff == 0 ? 0 : (j*(j+1.0)-l*(l+1.0)-m*(m+1.0))
  *(coeff/(R*a))*(pow(hbarc/pion_mass,2))
  *(exp((R-Rr)/a))/(pow((1.0+exp((R-Rr)/a)),2.0)));

  return Vspin;
}

double Gauss_potential::getValue(double R)const
{
  return (coeff == 0 ? 0 : -1.0*coeff*exp(-1.0*pow(R/r,2)));
}

double CouplingPotential::getValue(double r1, double r2, Particle p)const
{
  double Rr = r;//*pow(p.getM(),1.0/3.0);
  double R = (r1+r2)/2.0;
  double pot = (coeff == 0 ? 0 : beta_c//*r/a
    *(coeff*exp((R-Rr)/a))/pow((1.0+exp((R-Rr)/a)),2.0));
  if(beta_nl != 0 && coeff != 0) pot *= //exp(-1.0*(r1*r1 + r2*r2)/(beta_nl*beta_nl));
    exp(-1.0*pow((r1-r2)/beta_nl, 2.0));
  else if(r1 != r2) return 0;
  return pot;
}
