#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <complex>
#include <stdint.h>
#include <string.h>

using namespace std;




//min and max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) > (b)) ? (b) : (a))

const int NOTHING = 0;
const double precision = 1E-10;
const double INFINITE = 1E250;
const complex<double> NICHTS =complex<double > (1E-150,1E-150);


//const double two_amu_over_hbar_square = 0.0478450427;
const double two_amu_over_hbar_square = 0.0478450697;

const double Coulomb_constant = 1.4399651788;
//const double proton_mass = 1.00783;


const double neutron_mass = 1.00866;
const double proton_mass = 1.0072764;

const double diproton_mass = 2.01566;
const double electron_mass = 0.510998902;
//const double amu = 931.49432;
const double amu = 931.494061;


const double fine_struct_const = 0.00729735253;

//const double hbar_c = 197.327053;

const double hbar_c = 197.32697;

const double Compton_wavelength_electron = 386.1592642;

const double cos_theta_tab[4] = {M_SQRT1_2,-M_SQRT1_2,-M_SQRT1_2, M_SQRT1_2};
const double sin_theta_tab[4] = {M_SQRT1_2, M_SQRT1_2,-M_SQRT1_2,-M_SQRT1_2};

// two_amu_over_hbar_square : 2amu/(h/2pi)²; amu : atomic mass unit 
// cst_coulomb : the 1.44 MeV in the Coulomb potential. 
// proton_mass, neutron_mass : masses in amu 

enum particle_type {NEUTRON,PROTON,DIPROTON};
enum problem_type {E_SEARCH,VO_SEARCH,VSO_SEARCH,E_SCAT};
enum interaction_type {FULL_DELTA,CUT_DELTA,SDI,SKYRME,GOGNY,CM,P2_OVER_2M,FROM_FILE};

// Functions integrated in the Schrödinger equation.
enum potential_type {COULOMB,WS,PTG,SELF_CONSISTENT,HF_EFF_MASS,HF_NOT_EFF_MASS,NO_POT};


 double sqrt( int a)
{
  return sqrt(a*1.0);

}


complex<double> operator * ( complex<double> z, int a)
{
  double b=a*1.0;

  return complex<double> (real (z)*b,imag (z)*b);
}

complex<double> operator * ( int a, complex<double> z)
{
  double b=a*1.0;

  return complex<double> (real (z)*b,imag (z)*b);
}




//sign function.
inline char SIGN (const double a)
{
  if (a < 0)
    return -1;
  else 
    return 1;
}

// integer from double
inline int make_int (const double x_int)
{ 
  return static_cast<int> (rint (x_int));
}

// unsigned integer from double
inline unsigned int make_uns_int (const double x_int)
{
  return static_cast<unsigned int> (rint (x_int));
}


// unsigned char from double
inline unsigned char make_uns_char (const double x_int)
{
  return static_cast<unsigned char> (rint (x_int));
}


//-1 pow integer from double
inline int minus_one_pow(const double x_int)
{
  return ((make_int (x_int) % 2 == 0) ? (1) : (-1));
}

double real (const double x)
{
  return x;
}

double imag (const double x)
{
  return 0.0;
}


double inf_norm (const complex<double> z)
{
  return max (abs (real (z)),abs (imag (z)));
}
                                                                                                             
double inf_norm (const double x)
{
  return abs (x);
}



double gamma(double x)
{
    int i,k,m;
    double ga,gr,r,z;

    static double g[] = {
        1.0,
        0.5772156649015329,
       -0.6558780715202538,
       -0.420026350340952e-1,
        0.1665386113822915,
       -0.421977345555443e-1,
       -0.9621971527877e-2,
        0.7218943246663e-2,
       -0.11651675918591e-2,
       -0.2152416741149e-3,
        0.1280502823882e-3,
       -0.201348547807e-4,
       -0.12504934821e-5,
        0.1133027232e-5,
       -0.2056338417e-6,
        0.6116095e-8,
        0.50020075e-8,
       -0.11812746e-8,
        0.1043427e-9,
        0.77823e-11,
       -0.36968e-11,
        0.51e-12,
       -0.206e-13,
       -0.54e-14,
        0.14e-14};

    if (x > 171.0)
      {
	return 1e308;  
	
	cout<<"x too large"<<endl;
	exit(1);
      }


    if (x == (int)x) {
        if (x > 0.0) {
            ga = 1.0;               // use factorial
            for (i=2;i<x;i++) {
               ga *= i;
            }
         }
         else
            ga = 1e308;
     }
     else {
        if (fabs(x) > 1.0) {
            z = fabs(x);
            m = (int)z;
            r = 1.0;
            for (k=1;k<=m;k++) {
                r *= (z-k);
            }
            z -= m;
        }
        else
            z = x;
        gr = g[24];
        for (k=23;k>=0;k--) {
            gr = gr*z+g[k];
        }
        ga = 1.0/(gr*z);
        if (fabs(x) > 1.0) {
            ga *= r;
            if (x < 0.0) {
                ga = -M_PI/(x*ga*sin(M_PI*x));
            }
        }
    }
    return ga;
}
