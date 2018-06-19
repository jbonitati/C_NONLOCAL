#include <cmath>
#include "particle.h"
#include "constants.h"

/*Returns the Coulomb potential for particles at distance R*/
double Particle::coulomb_potential_to(double R, double rc, Particle targ)const
{
  //return 6*E2HC / R;
  
  double Rcoul = rc;//*(pow(getM(), 1.0/3.0) + pow(targ.getM(),1.0/3.0));
  double z1 = targ.getZ();
  double z2 = getZ();

  if (R<Rcoul)
  {
    return ((z1*z2*E2HC)/(2.0))*(3.0-pow(R,2.0)/(pow(Rcoul,2)));
  }
  else
  {
    return (z1*z2*E2HC)/R;
  }
}  
