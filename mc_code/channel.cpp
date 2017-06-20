#include <cmath>
#include "channel.h"
#include "particle.h"

const double hbarc = 197.3269718;
const double mass_unit=931.494;

/*Returns the potential due to the centrifugal term in the Hamiltonian*/
double Channel::central_potential(double R)
{
	double constant = (pow(hbarc,2)/(2.0*mass_unit*mu));

	return constant*((l*(l+1.0))/(R*R));
}

