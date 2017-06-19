#include <complex>
#include <cmath> //for pow and exp
#include <gsl/gsl_integration.h>
#include "potential_functions.h"

using namespace std;

const double hbarc = 197.3269718;
const double mass_unit=931.494;
const double c_constant=(pow(hbarc,2)/(2.0*mass_unit));
const double PI = 3.14159265359;
const std::complex<double> I(0., 1.);

/*Returns the Coulomb potential for particles with proton numbers z1,z2
 *  at distance R*/
double coulomb_potential(int z1,int z2,double R,int n, double rc)
{
    double Rcoul,Vcoul;
    Rcoul=rc*pow((n+z1),1.0/3.0);
  
    if (R<Rcoul)
	  {
		Vcoul=((z1*z2*1.43997)/(2.0*Rcoul))*(3.0-pow(R,2.0)/(pow(Rcoul,2)));
	  }
    else
      {
		Vcoul=(z1*z2*1.43997)/R;
      }
    return Vcoul;
}  

/*Returns the potential due to the centrifugal term in the Hamiltonian*/
double central_potential(int l,double R,double u)
{
	double constant;

	constant=c_constant/u;

	return constant*((l*(l+1.0))/(R*R));
}


/* Vd and Wd are real and imaginary surface potentials
 * V and Wv are real and imaginary volume potentials
 * Vsp and Wsp are real and imaginary spin orbit potentials */

double Vd_potential(double Vd, double R,int n,int z,double rvd,double avd)
{
    double rvvd;
    rvvd=rvd*pow(n+z,1.0/3.0);
    return (-4.0*Vd*exp((R-rvvd)/avd))/pow((1.0+exp((R-rvvd)/avd)),2.0) ;
}

double V_sp_potential(double Vso,double Rso,double aso,
	double r,double j,int l,int n,int z)
{
   double Vspin,Rrso;
   Rrso=Rso*pow(n+z,1.0/3.0);
   Vspin=(j*(j+1.0)-l*(l+1.0)-0.5*(0.5+1.0))
	*(-Vso/(aso*r))*(pow(hbarc/139.6,2))
	*(exp((r-Rrso)/aso))/(pow((1.0+exp((r-Rrso)/aso)),2.0));

   return Vspin;
}

double W_sp_potential(double Wso,double Rwso,double awso,
	double r,double j,int l, int n, int z)
{
   double Wspin,Rwwso;
   Rwwso=Rwso*pow(n+z,1.0/3.0);
   Wspin=(j*(j+1.0)-l*(l+1.0)-0.5*(0.5+1.0))
	*(-Wso/(awso*r))*(pow(hbarc/139.6,2))
	*(exp((r-Rwwso)/awso))/(pow((1.0+exp((r-Rwwso)/awso)),2.0));
 
   return Wspin;
}

double  Wd_potential(double Wd,double R,int n,int z,double rwd,double awd)
{
    double rwwd;
    rwwd=rwd*pow(n+z,1.0/3.0);
 
    return (-4.0*Wd*exp((R-rwwd)/awd))/pow((1.0+exp((R-rwwd)/awd)),2.0) ;
}


double  f(double R,int n,int z,double r,double a)
{   
    double Rr;
    Rr=r*pow(n+z,1.0/3.0);

    return 1.0/(1.0+exp((R-Rr)/a))  ;
}


double Wv_potential(double Wv,double R,int n,int z,double rwv,double awv)  
{
    return -Wv*f(R,n,z,rwv,awv);
}


double  V_potential(double Vv,double R,int n,int z,double rv,double av)
{
    return -Vv*f(R,n,z,rv,av);
}


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



/*Calculate the total Woods-Saxon potential*/
std::complex<double> local_potential(double R,int l,double j,int n,
	int z1,int z2,double Vv,double rv,double av,double Wv,
	double rwv,double awv,double Vd,double rvd,double avd,double Wd,
	double rwd,double awd,double Vso,double Rso,double aso,double Wso,
	double Rwso,double awso)
{   
    return V_potential(Vv,R,n,z1,rv,av)
		+1.0*I*Wv_potential(Wv,R,n,z1,rwv,awv)
		+Vd_potential(Vd,R,n,z1,rvd,avd)
		+1.0*I*Wd_potential(Wd,R,n,z1,rwd,awd)
		+V_sp_potential(Vso,Rso,aso,R,j,l,n,z1)
		+1.0*I*W_sp_potential(Wso,Rwso,awso,R,j,l,n,z1)
		+coulomb_potential(z1,z2,R,n, 1.25);
}


std::complex<double> non_local(double r_minus,double r_plus,int l,
	double jj,int n,int z1,int z2, double beta,double Vv1,double rv1,
	double av1,double Wv1,double rwv1,double awv1,double Vd1,double rvd1,
	double avd1,double Wd1,double rwd1,double awd1,double Vso1,
	double Rso1,double aso1,double Wso1,double Rwso1,double awso1)
{
    double r_average,zz;
  
    r_average=(r_plus+r_minus)/2.0;
     /*
    complex < double> c_n,c,total;
    zz=2.0*r_minus*r_plus/(beta*beta);
    c_n=2.0*pow(1.0*I,l)*zz;
    c=1.0/(2.0*pow(1.0*I,l));
    complex <double> V_potential_average;
    */
    V_potential_average=local_potential(r_average,l,jj,n,z1,z2,Vv1,
		rv1,av1,Wv1,rwv1,awv1,Vd1,rvd1,avd1,Wd1,rwd1,awd1,Vso1,Rso1,
		aso1,Wso1,Rwso1,awso1);
    /*
    if ( abs(zz)<=700)
    {  
		 complex <double> part_b,part_a;
     
         part_a=c*integration(l,zz);

         complex <double> kl;
  
         double exponent;
         kl=c_n*part_a;
       
         
         exponent=(pow(r_plus,2.0)+pow(r_minus,2.0))/(beta*beta);
         total=1.0/(beta*pow(PI,1.0/2.0))*(exp(-exponent))*kl ;
    }
    else
    {
         total=1.0/(beta*pow(PI,1.0/2.0))*exp(-pow((r_plus-r_minus)/beta,2.0));
    }
   */
   total = 1.0/(pow(beta,3.0)*pow(PI,3.0/2.0))*exp(-pow((r_plus-r_minus)/beta,2.0));
    return V_potential_average*total;
   // return 0;
}
