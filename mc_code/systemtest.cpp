//based on tutorial from http://neyasystems.com/an-engineers-guide-to-unit-testing-cmake-and-boost-unit-tests/
//Link to Boost
 #define BOOST_TEST_DYN_LINK

//Define our Module name (prints at testing)
 #define BOOST_TEST_MODULE "BaseClassModule"

//VERY IMPORTANT - include this last
#include <boost/test/unit_test.hpp>
#include "system.h"
#include <iostream>

using namespace std;

// ------------- Tests Follow --------------
//Name your test cases for what they test
BOOST_AUTO_TEST_CASE( constructor )
{
  double a = 10.0;
  double energy = 5.0;
  Particle projectile(1,0);
  Particle target(10,4);
  double rc = 1.0; //coulomb radius
  int c0 = 1; //entrance channel
  OpticalPotential op(V_Potential(0, 0, 0), V_Potential(0, 0, 0),
    D_Potential(0, 0, 0), D_Potential(0,0,0),
    SO_Potential(0, 0, 0), SO_Potential(0, 0, 0));
  NonLocalOpticalPotential nlop(V_Potential(0, 0, 0), 
    V_Potential(0, 0, 0),
    D_Potential(0, 0, 0), D_Potential(0,0,0),
    SO_Potential(0, 0, 0), SO_Potential(0, 0, 0), 0);
  int b_size = 10; //basis size
  double step = .1;
  double max = 15;
  matrix<CouplingPotential> coupling;
  std::vector<Channel> channels;
  double mu = 1.0;//reduced mass
  channels.push_back(Channel(
    0,//e
    0,//l
    0,//m
    0,//j
    mu, energy, a, target, projectile));
      
  //One channel case
  System system1(const double a, double energy,
    Particle projectile, Particle target, double rc,
    int c0, std::vector<Channel> c,
    OpticalPotential op, NonLocalOpticalPotential nlop, int b_size,
    double step, double max, matrix<CouplingPotential> coupling);
  
  //////////////////////////////////////////////////////////////////////
  //todo: add tests here, e.g:
  //BOOST_EQUALS(system1.Cmatrix, ((0, 1), (1, 0)));
  
  std::cout << "Testing!" << std::endl;
}
