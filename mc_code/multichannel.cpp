#include "system.h"

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>

#include <armadillo>

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstring>
#include <string>
#include <iomanip>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/lexical_cast.hpp>


const double hbarc = 197.3269718;
const double mass_unit=931.494;
const double c_constant=(pow(hbarc,2)/(2.0*mass_unit));
const int MAX_PTS = 6000;
const double E2HC = 0.00729927;
const double PI = 3.14159265359;
const std::complex<double> I(0., 1.);

//namespace is superior to static
namespace { 
  System * mySystem;
  std::string filename;
}

void loadSystem(){
  
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini("config.ini", pt);
  
  std::string filename = pt.get<std::string>("Settings.output_file");
  int num_channels = pt.get<int>("Settings.num_channels");
  int c0 = pt.get<int>("Settings.entrance_channel");
  
  std::vector<Channel*> channels;
  channels.reserve(num_channels);
  
  try{
    for(int i = 1; i <= num_channels; i++){
    channels.push_back(new Channel(
      pt.get<double>("Channel" + boost::lexical_cast<std::string>(i) + ".Energy"), 
      pt.get<double>("Channel" + boost::lexical_cast<std::string>(i) + ".Angular_momentum"), 
      pt.get<double>("Channel" + boost::lexical_cast<std::string>(i) + ".Total_angular_momentum"),
      pt.get<double>("Channel" + boost::lexical_cast<std::string>(i) + ".Spin")));
    }
  }catch(...){
    std::cerr << "\nError reading channel values. Make sure all channels are specified in config" << std::endl << std::endl;
    throw;
  }
  int NN=pt.get<int>("Numerical.Basis_size") ;

  double m1=pt.get<double>("Numerical.Projectile_mass_number") ;
  double m2=pt.get<double>("Numerical.Target_mass_number") ;
  double a_size=pt.get<double>("Numerical.Channel_radius") ;
  int z1=pt.get<int>("Numerical.Target_proton_number");
  int z2=pt.get<int>("Numerical.Projectile_proton_number");
  double Nr=pt.get<double>("Numerical.Step_size");
  double E = pt.get<double>("Numerical.Projectile_energy");
  double R_max=pt.get<double>("Numerical.R_max");

  double Vv=pt.get<double>("local.Vv") ;
  double rv=pt.get<double>("local.rv") ;
  double av=pt.get<double>("local.av") ;
  double Wv=pt.get<double>("local.Wv") ;
  double rwv=pt.get<double>("local.rwv") ;
  double awv=pt.get<double>("local.awv");
  double Vd=pt.get<double>("local.Vd");
  double rvd=pt.get<double>("local.rvd");
  double avd=pt.get<double>("local.avd");
  double Wd=pt.get<double>("local.Wd") ;
  double rwd=pt.get<double>("local.rwd");
  double awd=pt.get<double>("local.awd");
  double Vso=pt.get<double>("local.Vso");
  double Rso=pt.get<double>("local.Rso");
  double aso=pt.get<double>("local.aso");
  double Wso=pt.get<double>("local.Wso") ;
  double Rwso=pt.get<double>("local.Rwso");
  double awso=pt.get<double>("local.awso");
  
  OpticalPotential op(V_Potential(Vv, rv, av), V_Potential(Wv, rwv, awv),
    D_Potential(Vd, rvd, avd), D_Potential(Wd,rwd,awd),
    SO_Potential(Vso, Rso, aso), SO_Potential(Wso, Rwso, awso));

  double Vv1=pt.get<double>("Non_local.Vv") ;
  double rv1=pt.get<double>("Non_local.rv") ;
  double av1=pt.get<double>("Non_local.av") ;
  double Wv1=pt.get<double>("Non_local.Wv") ;
  double rwv1=pt.get<double>("Non_local.rwv") ;
  double awv1=pt.get<double>("Non_local.awv");
  double Vd1=pt.get<double>("Non_local.Vd");
  double rvd1=pt.get<double>("Non_local.rvd");
  double avd1=pt.get<double>("Non_local.avd");
  double Wd1=pt.get<double>("Non_local.Wd") ;
  double rwd1=pt.get<double>("Non_local.rwd");
  double awd1=pt.get<double>("Non_local.awd");
  double Vso1=pt.get<double>("Non_local.Vso");
  double Rso1=pt.get<double>("Non_local.Rso");
  double aso1=pt.get<double>("Non_local.aso");
  double Wso1=pt.get<double>("Non_local.Wso") ;
  double Rwso1=pt.get<double>("Non_local.Rwso");
  double awso1=pt.get<double>("Non_local.awso");

  double beta=pt.get<double>("Non_local.beta");
  
  double coupling = pt.get<double>("coupling.beta");
  
  NonLocalOpticalPotential nlop(V_Potential(Vv1, rv1, av1), 
    V_Potential(Wv1, rwv1, awv1),
    D_Potential(Vd1, rvd1, avd1), D_Potential(Wd1,rwd1,awd1),
    SO_Potential(Vso1, Rso1, aso1), SO_Potential(Wso1, Rwso1, awso1), beta);
    
  mySystem = new System(a_size, E, m1, m2, z1, z2, 
    num_channels, c0, channels, op, nlop, NN, Nr, R_max, coupling);
}

int main()
{   
  loadSystem();
  
  std::ofstream outfile(("output/" + filename).c_str());
  
  mySystem->waveFunction(outfile);
  
  outfile.close();
  
  return 0;
}
