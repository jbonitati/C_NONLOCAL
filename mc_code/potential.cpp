#include <cmath>
#include "potential.h"


const double hbarc = 197.3269718;

double V_Potential::getValue(double R, Particle p)const{
  double Rr = r*pow(p.getM(),1.0/3.0);
  return coeff*-1.0/(1.0+exp((R-Rr)/a));
}

double D_Potential::getValue(double R, Particle p)const{
  double Rr = r*pow(p.getM(),1.0/3.0);
  return (-4.0*coeff*exp((R-Rr)/a))/pow((1.0+exp((R-Rr)/a)),2.0);
}

double SO_Potential::getValue(double R, Particle p, Channel * c)const{
  double Rr = r*pow(p.getM(),1.0/3.0);
  double l = c->getL();
  double j = c->getJ();
  double m = c->getSpin();
  
  double Vspin=(j*(j+1.0)-l*(l+1.0)-m*(m+1.0))
  *(-coeff/(a*R))*(pow(hbarc/139.6,2))
  *(exp((R-Rr)/a))/(pow((1.0+exp((R-Rr)/a)),2.0));

  return Vspin;
}

double Coupling_Potential::getValue(double r1, double r2, Particle p)const{
  double Rr = r*pow(p.getM(),1.0/3.0);
  double R = (r1+r2)/2.0;
  return beta_c*(-1.0*coeff*exp((R-Rr)/a))/pow((1.0+exp((R-Rr)/a)),2.0)
    *exp(-1.0*pow((r1-r2)/beta_nl, 2.0));
}
