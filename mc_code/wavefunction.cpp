#include "wavefunction.h"

//change this to print real/imaginary parts or both
std::ostream& operator<<(std::ostream &os, const WaveFunction &wf)
{
  //prints "[radius] [Channel 1 wave function] [Channel 2 wave function] ..."
  os << wf.r;
  for(std::vector<cx_double>::const_iterator it = wf.values.begin(); it != wf.values.end(); it++)
  {
    os << " " 
    << std::real((*it))
    << " "
    << std::imag((*it));
    //<< " "
    //<< std::abs((*it));
  }
  os << std::endl;
  return os;
}
