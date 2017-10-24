//based on tutorial from http://neyasystems.com/an-engineers-guide-to-unit-testing-cmake-and-boost-unit-tests/
//Link to Boost
 #define BOOST_TEST_DYN_LINK

//Define our Module name (prints at testing)
 #define BOOST_TEST_MODULE "BaseClassModule"

//VERY IMPORTANT - include this last
#include <boost/test/unit_test.hpp>
#include "system.h"
#include <iostream>

// ------------- Tests Follow --------------
//Name your test cases for what they test
BOOST_AUTO_TEST_CASE( constructor )
{
  /*System system1(const double a_size, double e,
      Particle projectile, Particle target, double rc,
      int c0, std::vector<Channel> channels_,
      OpticalPotential op, NonLocalOpticalPotential nlop, int basis_size,
      double step, double max,matrix<CouplingPotential> coupling);
  */
  
  std::cout << "Testing!" << std::endl;
}
