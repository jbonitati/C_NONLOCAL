#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <fenv.h>   // enable floating point trap
#include <fstream>
#define MAXN 1000
#define MAXX0 1000
#include <complex>
#include <cmath>
#include <vector>
#include <math.h>       /* pow */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
using namespace std;
static void __attribute__((constructor)) trapfpe () {
 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}  // program will stop if an invalid number is generated


double Legendre(int n, double t) // return P_{n}(t)
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

double f  (double x, void * params) 
{ 
  //const complex<double> i(0.0,1.0);
  //complex<double> c;   
  // c=1.0/(2.0*pow(1.0*i,l));
  double mu = *(double *) params;
  int l=*(int *) params;
 // cout<<"ppp";
 // cout<<exp(1.0*0.3)*Legendre(1, 0.3);
  return  exp(mu*x)*Legendre(l, x);
}


// helloacm.com
double integral(double(*f)(double x), double a, double b, int n) {
    double step = (b - a) / n;  // width of each small rectangle
    double area = 0.0;  // signed area
    for (int i = 0; i < n; i ++) {
        area += f(a + (i + 0.5) * step) * step; // sum up each small rectangle
    }
    return area;
}

double f_f(double x,int l, double mu)
{
 return  exp(mu*x)*Legendre(l, x);
}

double simpsons(int n, double a, double b,int l, double mu)
{
int i;
double dx, x, sum;

dx = (b - a) / n;
sum = f_f(a,l,mu) + f_f(b,l,mu);
for (i = 1; i < n; i++) {
x = a + dx * i;
sum += 2.0 * (1.0 + i%2) * f_f(x,l,mu);
}
sum *= dx/3.0;
return sum;
}

int
main (void)
{ 
  double result, error;
  double expected = -2.0;
  double mu = 4.0;
  int l=3;
  
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (2000);

  gsl_function F;
  F.function = &f;
  F.params = &mu;
  F.params = &l;
 // integral(double(*f)(double x), double a, double b, int n) 
  gsl_integration_qags (&F, -1.0, 1.0, 1e-7, 1e-7, 2000,
                        w, &result, &error); 
//  cout << integral(*f, -1, 1, 10);
  cout<<simpsons(12000, -1.0, 1.0,3, 4);
  
// cout<<Legendre(2, -0.3); 
  cout<<" ";
  cout<<result;
  cout<<" ";
  cout<<error;
  cout<<" ";
  gsl_integration_workspace_free (w);

  return 0;
}
