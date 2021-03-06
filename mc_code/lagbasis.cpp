#include <gsl/gsl_sf.h>
#include <cmath>
#include <iostream>
#include "lagbasis.h"

//Weichuan's code for computing the first n weights and abcissae for the 
//Gauss Legendre quadrature rule used in calculating lagrange functions
void LagBasis::gauleg(const double x1, const double x2, int n)
{
	const double EPS=1.0e-10;
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=0;i<m;i++) 
    {
		z=cos(3.141592654*(i+0.75)/(n+0.5));
		do 
        {
			p1=1.0;
			p2=0.0;
			for (j=0;j<n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
	
        xi[i]=(xm-xl*z+1.0)/2.0;
		xi[n-1-i]=(xm+xl*z+1.0)/2.0;
    //removed a factor of 2 from weights
    //in previous program, this was factored out when calculating potential
		weights[i]=xl/((1.0-z*z)*pp*pp);
		weights[n-1-i]=weights[i];
	}
}

//returns the value phi_i(r) where phi_i is the i^th lagrange function
//must have |2*r/a - 1| <= 1
double LagBasis::phi(int i, double r)const
{
  int N = size;
  if (r==a*xi[i-1])
    return pow(-1.0,N+i)*sqrt(a*xi[i-1]*(1.0-xi[i-1]));
  else 
    return pow(-1.0,N+i)*r/(a*xi[i-1])*sqrt(a*xi[i-1]*(1.0-xi[i-1]))
      *gsl_sf_legendre_Pl(N,2.0*r/a-1.0)/(r-a*xi[i-1]);
}

//returns a row vector conatining the basis functions evaluated at r
arma::rowvec LagBasis::get_phi_r(double r)const
{
  arma::rowvec phi_r(size);
  for(int i = 0; i < size; i++)
  {
    phi_r(i) = phi(i+1,r);
  }
  return phi_r;
}

//Lagrange Basis constructor
//precomputes all necessary values of the first n basis functions
LagBasis::LagBasis(int n, double a_): size(n), a(a_), phi_a(n)
{
  weights = new double[n];
  xi = new double[n];
  std::cout << "Precomputing basis functions..." << std::endl;
  gauleg(-1,1,n);
  for(int i = 0; i < n; i++)
  {
    phi_a(i) = phi(i+1,a);
    //std::cout << i << " " << phi(i+1,a) << std::endl;
  }
}

