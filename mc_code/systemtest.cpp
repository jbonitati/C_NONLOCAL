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
using namespace boost::unit_test;

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

// ------------- Tests Follow --------------
//Name your test cases for what they test
BOOST_AUTO_TEST_CASE( constructor, * utf::tolerance(0.0001) )
{
  std::cout << "Testing constructor!" << std::endl;

  
  double a = 10.0;
  double energy = 5.0;
  double m1 = 1;
  double m2 = 10;
  double z1 = 0;
  double z2 = 4;
  Particle projectile(m1,z1);
  Particle target(m2,z2);
  double rc = 0.5; //coulomb radius
  int c0 = 1; //entrance channel
  OpticalPotential op(V_Potential(53.28, 2.483, 0.5), V_Potential(0, 0, 0),
    D_Potential(0, 0, 0), D_Potential(0,0,0),
    SO_Potential(0, 0, 0), SO_Potential(0, 0, 0));
  NonLocalOpticalPotential nlop(V_Potential(0, 0, 0), 
    V_Potential(0, 0, 0),
    D_Potential(0, 0, 0), D_Potential(0,0,0),
    SO_Potential(0, 0, 0), SO_Potential(0, 0, 0), 0);
  int b_size = 10; //basis size
  double step = .1;
  double max = 20;
  matrix<CouplingPotential> coupling;
  std::vector<Channel> channels;
  double mu = (m1*m2)/(m1+m2);//reduced mass
  channels.push_back(Channel(
    0,//e
    0,//l
    0,//m
    0,//j
    mu, energy, a, target, projectile));
  cout << channels.size()<< endl;    
  
  //One channel case
  cout<< "Testing one channel case" << endl;
  System system1(a, energy, projectile,target, rc,
    c0, channels, op, nlop, b_size, step, max, coupling);
  //////////////////////////////////////////////////////////////////////
  //todo: add tests here, e.g:
  //BOOST_EQUALS(system1.Cmatrix, ((0, 1), (1, 0)));
  cout << system1.internalWaveFunction(0,10) << endl;
  //BOOST_TEST(std::real(system1.getRmatrix()(0,0)) == 0.8418078);
}


/*
  double a = 10.0;
  double energy = 5.0;
  double m1 = 1;
  double m2 = 10;
  double z1 = 0;
  double z2 = 4;
  Particle projectile(m1,z1);
  Particle target(m2,z2);
  double rc = 0.5; //coulomb radius
  int c0 = 1; //entrance channel
  OpticalPotential op(V_Potential(53.28, 2.483, 0.5), V_Potential(0, 0, 0),
    D_Potential(0, 0, 0), D_Potential(0,0,0),
    SO_Potential(0, 0, 0), SO_Potential(0, 0, 0));
  NonLocalOpticalPotential nlop(V_Potential(0, 0, 0), 
    V_Potential(0, 0, 0),
    D_Potential(0, 0, 0), D_Potential(0,0,0),
    SO_Potential(0, 0, 0), SO_Potential(0, 0, 0), 0);
  int b_size = 10; //basis size
  double step = .1;
  double max = 20;
  matrix<CouplingPotential> coupling;
  std::vector<Channel> channels;
  double mu = (m1*m2)/(m1+m2);//reduced mass
  channels.push_back(Channel(
    0,//e
    0,//l
    0,//m
    0,//j
    mu, energy, a, target, projectile));

  //One channel case
  cout<< "Testing one channel case" << endl;
  System system1(a, energy, projectile,target, rc,
    c0, channels, op, nlop, b_size, step, max, coupling);
*/
