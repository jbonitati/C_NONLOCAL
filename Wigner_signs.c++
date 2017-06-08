//$Id: Wigner_signs.c++,v 1.6 2002/10/01 16:39:37 lmichel Exp $

// Factorial n! for a double x=n which is in fact an integer.
double factorial (const double x)
{
  const int n = make_int (x);
  double fact = 1.0;

  for (int i = 1 ; i <= n ; i++)
    fact *= i;

  return fact;
}

// Factorial n!! for a double x=n which is in fact an integer.
double double_factorial (const double x)
{
  const int n = make_int (x);      
  double double_fact = 1.0;

  for (int i = n ; i > 0 ; i -= 2)
    double_fact *= i;

  return double_fact;
}


// Test to know if the three j's verify the triangle relation.

bool is_it_triangle (const double j1,const double j2,const double j3)
{
  return ((j3 >= abs (j1 - j2)) && (j3 <= j1 + j2));
}

// Test to know if a double is a integer.

bool is_it_integer (double x)
{
  return (abs (x - rint (x)) < precision);
}

bool triang (const double j1,const double j2,const double j3)
{
  return ((j3 >= abs (j1 - j2)) && (j3 <= j1 + j2) && is_it_integer(j1+j2+j3));
}


// Square root of the delta appearing in the Wigner signs.

double sqrt_delta (const double j1,const double j2,const double j3)
{
  double log_delta_term;

  log_delta_term = 
     log (factorial (j1 + j2 - j3))
     + log (factorial (j2 + j3 - j1))
     + log (factorial (j3 + j1 - j2))
     - log (factorial (j1 + j2 + j3 + 1.0));

  log_delta_term *= 0.5;

  return (exp (log_delta_term));
}


// Wigner 3j.
// The Racah sum is used to calculate it : formula (22) p.910 in Messiah 2.

double Wigner_3j (const double j1,const double j2,const double j3,
		  const double m1,const double m2,const double m3)
{
  if (rint (m1 + m2 + m3) != 0.0) return 0.0;
  if (is_it_triangle (j1,j2,j3) == false) return 0.0;
  if (is_it_integer (j1 - m1) == false) return 0.0;
  if (is_it_integer (j2 - m2) == false) return 0.0;
  if (is_it_integer (j3 - m3) == false) return 0.0;
  
  const int tmin0 = make_int (max (0.0,j2 - j3 - m1));
  const int tmin = make_int (max (tmin0,j1 - j3 + m2));

  const int tmax0 = make_int (min (j1 + j2 - j3,j1 - m1));
  const int tmax = make_int (min (tmax0,j2 + m2));

  double Wigner_sign = 0.0;
  for (int t = tmin ; t <= tmax ; t++)
    Wigner_sign += minus_one_pow (t)
      /factorial (t)
      /factorial (j3 - j2 + m1 + t)
      /factorial (j3 - j1 - m2 + t)
      /factorial (j1 + j2 - j3 - t)
      /factorial (j1 - m1 - t)
      /factorial (j2 + m2 - t);
    
  Wigner_sign *= minus_one_pow (j1 - j2 - m3)
    *sqrt (factorial (j1 + m1)*factorial (j1 - m1)
	  *factorial (j2 + m2)*factorial (j2 - m2)
	  *factorial (j3 + m3)*factorial (j3 - m3))
    *sqrt_delta (j1,j2,j3);

  if (abs (Wigner_sign) < precision) Wigner_sign = 0.0;
  
  return Wigner_sign;
}

// Clebsch - Gordan
// It is calculated with the corresponding 3j.

double Clebsch_Gordan (const double j1,const double m1,
		       const double j2,const double m2,
		       const double J ,const double M)
{
  return sqrt (2.0*J + 1.0)*minus_one_pow (j1 - j2 + M)*Wigner_3j (j1,j2,J,m1,m2,-M);
}


// Wigner 6j.
// The Racah sum is used to calculate it : formula (36) p.915 in Messiah 2.

double Wigner_6j (const double j1,const double j2,const double j3,
		  const double J1,const double J2,const double J3)
{
  if (is_it_triangle (j1,j2,j3) == false) return 0.0;
  if (is_it_triangle (j1,J2,J3) == false) return 0.0;
  if (is_it_triangle (J1,j2,J3) == false) return 0.0; 
  if (is_it_triangle (J1,J2,j3) == false) return 0.0;
  if  (is_it_integer (j1 + j2 + j3) == false) return 0.0;
  if  (is_it_integer (j1 + J2 + J3) == false) return 0.0;
  if  (is_it_integer (J1 + j2 + J3) == false) return 0.0;
  if  (is_it_integer (J1 + J2 + j3) == false) return 0.0;

  const int tmin0 = make_int (max (j1 + j2 + j3,j1 + J2 + J3));
  const int tmin1 = make_int (max (tmin0,J1 + j2 + J3));
  const int tmin = make_int (max (tmin1,J1 + J2 + j3));

  const int tmax0 = make_int (min (j1 + j2 + J1 + J2,j2 + j3 + J2 + J3));
  const int tmax = make_int (min (tmax0,j3 + j1 + J3 + J1));

  double Wigner_sign = 0.0;
  for (int t = tmin ; t <= tmax ; t++)  
    Wigner_sign += minus_one_pow (t)
      *factorial (t + 1)
      /factorial (t - j1 - j2 - j3)
      /factorial (t - j1 - J2 - J3)
      /factorial (t - J1 - j2 - J3)
      /factorial (t - J1 - J2 - j3)
      /factorial (j1 + j2 + J1 + J2 - t)
      /factorial (j2 + j3 + J2 + J3 - t)
      /factorial (j3 + j1 + J3 + J1 - t);	                     
    
  Wigner_sign *= sqrt_delta (j1,j2,j3)*sqrt_delta (j1,J2,J3)
                *sqrt_delta (J1,j2,J3)*sqrt_delta (J1,J2,j3);
  
  if (abs (Wigner_sign) < precision) Wigner_sign = 0.0;
  
  return Wigner_sign;
}




// Wigner 9j
// It is calculated with 6j : formula (41) p.917 in Messiah 2.

double Wigner_9j (const double  j1,const double  j2,const double J12,
		  const double  j3,const double  j4,const double J34,
		  const double J13,const double J24,const double J)
{
  const double gmin0 = max(abs (j1 - J),abs (j3 - J24));
  const double gmin = max(gmin0,abs (j2 - J34));

  const double gmax0 = min(j1 + J,j3 + J24);
  const double gmax = min(gmax0,j2 + J34);

  double Wigner_sign = 0.0;
  for (double g = gmin ; rint (g - gmax) <= 0.0 ; g += 1.0)
    Wigner_sign += minus_one_pow (2.0*g)*(2.0*g + 1.0)
	          *Wigner_6j (j1,j2,J12,J34,J,g)
 	          *Wigner_6j (j3,j4,J34,j2,g,J24)
	          *Wigner_6j (J13,J24,J,g,j1,j3);

  if (abs (Wigner_sign) < precision) Wigner_sign = 0.0;

  return Wigner_sign;
}







