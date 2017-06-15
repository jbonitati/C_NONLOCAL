#include "potential_functions.h"

#include <iostream>
#include <cstdlib>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>

#include <armadillo>

#include <cstdio>
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <iomanip>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/lexical_cast.hpp>


const double hbarc = 197.3269718;
const double mass_unit=931.494;
const double c_constant=(pow(hbarc,2)/(2.0*mass_unit));
const int MAX_PTS = 6000;
const double E2HC = 0.00729927;
const double PI = boost::math::constants::pi<boost::multiprecision::cpp_dec_float_50>();
const std::complex<double> I(0., 1.);


/*Returns the value of the Legendre polynomial P_n at t*/
double Legendre(int n, double t)
{
	 int k;
	 double Pk_1,Pk_2,Pk; // P_{k-1}(x), P_{k-2}(x), P_k(x)

	 Pk_2 = 0.0;
	 Pk_1 = 1.0;
	 Pk = 1.0;
	 
	 for(k=1;k<=n;k++)
	 {
	  Pk = (2.0*k-1.0)/k*t*Pk_1 - (k-1.0)/k*Pk_2; 
	  Pk_2 = Pk_1;
	  Pk_1 = Pk;
	 }

	 return Pk;
}

struct my_f_params {int a; double b;};
double f1 (double x, void * p) 
{ //Derivative part of the integral
        struct my_f_params * params = (struct my_f_params *)p;
        int l = (params->a);
        double mu = (params->b);
        
        return exp(x*mu)*Legendre(l,x);
}
       
double integration(int l, double mu)
{
        std::complex<double> c;
        gsl_integration_workspace *work_ptr 
			= gsl_integration_workspace_alloc (1000);

       double lower_limit = -1.0;	/* lower limit a */
       double upper_limit = 1.0;	/* upper limit b */
       double abs_error = 1.0e-8;	/* to avoid round-off problems */
       double rel_error = 1.0e-8;	/* the result will usually be much better */
       double result;		/* the result from the integration */
       double error;			/* the estimated error from the integration */
    
        struct my_f_params alpha;
        gsl_function F1;       
        F1.function = &f1; 
        F1.params = &alpha;
        //cout<<"a"<<alpha.a<<endl;
        //c=1.0/(2.0*pow(1.0*i,alpha.a));
    
        alpha.a = l;
        alpha.b = mu;
        gsl_integration_qags (&F1, lower_limit, upper_limit,
			abs_error, rel_error, 1000, work_ptr, &result,
			&error);
        return result;
}





int main()
{   
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini("config.ini", pt);
  
  std::string out_name = pt.get<std::string>("Settings.output_file");
  int num_channels = pt.get<int>("Settings.num_channels");
  
  Channel* channels = new Channel[num_channels];
  
  for(int i = 1; i <= num_channels; i++){
    channels[i] = new Channel(
      pt.get<double>("Channel" + boost::lexical_cast<std::string>(i) + ".Energy"), 
      pt.get<double>("Channel" + boost::lexical_cast<std::string>(i) + ".Angular_momentum"), 
      pt.get<double>("Channel" + boost::lexical_cast<std::string>(i) + ".Total_angular_momentum")) 
  }
  
  int NN=pt.get<int>("Numerical.Basis_size") ;
  double m1=pt.get<double>("Numerical.Projectile_mass_number") ;
  double m2=pt.get<double>("Numerical.Target_mass_number") ;
  double a_size=pt.get<double>("Numerical.Channel_radius") ;
  int z1=pt.get<int>("Numerical.Target_proton_number");
  int z2=pt.get<int>("Numerical.Projectile_proton_number");
  int B=pt.get<int>("Numerical.BlochConstant") ;
  double Nr=pt.get<double>("Numerical.Step_size");
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

  double Vv1=pt.get<double>("Non_local.Vv1") ;
  double rv1=pt.get<double>("Non_local.rv1") ;
  double av1=pt.get<double>("Non_local.av1") ;
  double Wv1=pt.get<double>("Non_local.Wv1") ;
  double rwv1=pt.get<double>("Non_local.rwv1") ;
  double awv1=pt.get<double>("Non_local.awv1");
  double Vd1=pt.get<double>("Non_local.Vd1");
  double rvd1=pt.get<double>("Non_local.rvd1");
  double avd1=pt.get<double>("Non_local.avd1");
  double Wd1=pt.get<double>("Non_local.Wd1") ;
  double rwd1=pt.get<double>("Non_local.rwd1");
  double awd1=pt.get<double>("Non_local.awd1");
  double Vso1=pt.get<double>("Non_local.Vso1");
  double Rso1=pt.get<double>("Non_local.Rso1");
  double aso1=pt.get<double>("Non_local.aso1");
  double Wso1=pt.get<double>("Non_local.Wso1") ;
  double Rwso1=pt.get<double>("Non_local.Rwso1");
  double awso1=pt.get<double>("Non_local.awso1");

  double beta=pt.get<double>("Non_local.beta");
  
  Collision system = new Collision(a_size, B,m1, m2, z1, z2, channels);
  
  //calculate C matrix
  
  //invert C matrix
  
  //calculate R matrix
  
  //calculate U matrix
  
  //calculate wave function
  
  return 0;
}
