#include "particle.h"
#include <cmath>

/*Returns the Coulomb potential for particles at distance R*/
double Particle::coulomb_potential_to(double R, Particle targ)const
{
  double rc = 1.25; //don't know how to calculate this
  double Rcoul = rc*pow((targ.getM()),1.0/3.0);
  double z1 = targ.getZ();
  double z2 = getZ();

  if (R<Rcoul)
  {
    return ((z1*z2*1.43997)/(2.0))*(3.0-pow(R,2.0)/(pow(Rcoul,2)));
  }
  else
  {
    return (z1*z2*1.43997)/R;
  }
}  
