#include <iostream>
#include <iomanip>
#include <fstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include <gsl/gsl_poly.h>

using namespace std;


double kroneck(int i, int j)
{
  double a=0.0;

  if (i==j)
    {
      a=1.0;
    }

  return a;

}





int main()
{
  cout << "Hello World";        // prints Hello World

  return 0;
}
