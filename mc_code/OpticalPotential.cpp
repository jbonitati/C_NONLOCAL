#include <armadillo>
#include <cmath>
#include "OpticalPotential.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>

const arma::cx_double I(0., 1.);
const double PI = 3.14159265359;

arma::cx_double OpticalPotential::Volume_potential
(double radius, Particle p)const{
  return Volume_Re.getValue(radius, p) 
    + I*Volume_Im.getValue(radius, p);
}

arma::cx_double OpticalPotential::Surface_potential
(double radius, Particle p)const{
  return Surface_Re.getValue(radius, p) 
    + I*Surface_Im.getValue(radius, p);
}

arma::cx_double OpticalPotential::Spin_orbit_potential
(double radius, Particle p, Channel * c)const{
  return SpinOrbit_Re.getValue(radius, p, c)
    + I*SpinOrbit_Im.getValue(radius, p, c);
}

arma::cx_double OpticalPotential::totalPotential
(double radius, Particle p, Channel * c) const{
  return Volume_potential(radius, p) 
    + Surface_potential(radius, p) 
    + Spin_orbit_potential(radius, p, c);
}

struct my_f_params {int a; double b;};
 
double f1 (double x, void * p) 
{
  struct my_f_params * params = (struct my_f_params *)p;
  int l = (params->a);
  double mu = (params->b);

  return exp(x*mu)*gsl_sf_legendre_Pl(l,x);
}

double integration(int l, double mu)
{
  arma::cx_double c;
  gsl_integration_workspace *work_ptr 
    = gsl_integration_workspace_alloc (10000);

  double lower_limit = -1.0;	/* lower limit a */
  double upper_limit = 1.0;	/* upper limit b */
  double abs_error = 1.0e-7;	/* to avoid round-off problems */
  double rel_error = 1.0e-7;	/* the result will usually be much better */
  double result;		/* the result from the integration */
  double error;			/* the estimated error from the integration */

  struct my_f_params alpha;
  gsl_function F1;            
  F1.function = &f1; 
  F1.params = &alpha;
  //c=1.0/(2.0*pow(1.0*i,alpha.a));

  alpha.a = l;
  //cout<<"a"<<alpha.a<<endl;
  alpha.b = mu;
  gsl_integration_qags (&F1, lower_limit, upper_limit,
    abs_error, rel_error, 10000, work_ptr, &result,
    &error);
  return result;
}
arma::cx_double NonLocalOpticalPotential::totalPotential
(double r1, double r2, Particle p, Channel  * c) const{
  double ravg = (r1+r2)/2.0;
  int l = c->getL();
  
  arma::cx_double local_part 
    = OpticalPotential::totalPotential(ravg, p, c);
    
  double z =2.0*r1*r2/(beta*beta);
  arma::cx_double nonlocal_part;
  //for large z, we use an approximation
  if(abs(z) >= 500){
    nonlocal_part = 1.0/(pow(PI,0.5)*beta) * exp(-1.0*pow((r1-r2)/beta,2.0));
  }else{
    arma::cx_double approx_bessel,
      c_n=2.0*pow(1.0*I,l)*z,
      c_=1.0/(2.0*pow(1.0*I,l));

    approx_bessel =c_*integration(l,z);

    arma::cx_double kl;

    kl=c_n*approx_bessel;
    
    double exponent=(pow(r1,2.0)+pow(r2,2.0))/(beta*beta);
    nonlocal_part= kl * 1.0/(beta*pow(PI,1.0/2.0))*(exp(-1*exponent));
  }
    
  return local_part * nonlocal_part;
}
