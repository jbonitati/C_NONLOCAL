// ----------------------------------------------------------------------------
// 
// ----------------------------------------------------------------------------
// $Id: poly_splines.c++,v 1.6 2002/12/22 08:50:55 lmichel Exp $
// ----------------------------------------------------------------------------







// Calculation of a polynom with the Hörner scheme.
// ------------------------------------------------
// Class array version
//
// Variables:
// ----------
// degree : degree of the polynom.
// coeff : coefficients of the polynom : P(X) = \sum^{k=0}^{degree} coeff(k) X^k.
// x : variable of the polynom.
// Px : value of the polynom in x.

template <typename SCALAR_TYPE>
SCALAR_TYPE poly_eval (const unsigned int degree,
		       const class array<SCALAR_TYPE> &coeff,
		       const SCALAR_TYPE x)
{
  SCALAR_TYPE Px = coeff(degree);

  for (unsigned int i = degree ; i > 0 ; i--)
    Px = coeff(i-1) + x*Px;
  
  return Px;
}





// Calculation of a polynom with the Hörner scheme.
// ------------------------------------------------
// Standard table version
// 
// Variables:
// ----------
// degree : degree of the polynom.
// coeff : coefficients of the polynom : P(X) = \sum^{k=0}^{degree} coeff(k) X^k.
// x : variable of the polynom.
// Px : value of the polynom in x.

template <typename SCALAR_TYPE>
SCALAR_TYPE poly_eval (const unsigned int degree,
		       const SCALAR_TYPE coeff[],
		       const SCALAR_TYPE x)
{
  SCALAR_TYPE Px = coeff[degree];

  for (unsigned int i = degree ; i > 0 ; i--)
    Px = coeff[i-1] + x*Px;
  
  return Px;
}




// Calculation of a rational fraction with a polynom with the Hörner scheme.
// -------------------------------------------------------------------------
// The rational fonction is F = P/Q, with P and Q polynoms.
//
// Variables:
// ----------
// degree_num : degree of the polynom of the numerator.
// degree_den : degree of the polynom of the denominator.
// coeff_num : coefficients of the polynom of the numerator : P(X) = \sum^{k=0}^{degree_num} coeff_num(k) X^k.
// coeff_den : coefficients of the polynom of the denominator : Q(X) = \sum^{k=0}^{degree_den} coeff_den(k) X^k.
// x : variable of the polynom.
// Px : value of the polynom of the numerator in x.
// Qx : value of the polynom of the denominator in x.
// Fx : value of the rational fraction in x.

template <typename SCALAR_TYPE>
SCALAR_TYPE frac_eval (const unsigned int degree_num,
		       const class array<SCALAR_TYPE> &coeff_num,
		       const unsigned int degree_den,
		       const class array<SCALAR_TYPE> &coeff_den,
		       const SCALAR_TYPE x)
{
  const SCALAR_TYPE Px = poly_eval (degree_num,coeff_num,x);
  const SCALAR_TYPE Qx = poly_eval (degree_den,coeff_den,x);
  const SCALAR_TYPE Fx = Px/Qx;

  return Fx;
}






// Function He_smooth giving a smoothed Heaviside function.
// -----------------------------------------------------
// One has this function as it is better to have a function going smoothly from 0 to 1
// for integrations of integrals and differentials equations.
// The function is, with x0=-1 and x1=1, written He_smooth(x) :
// He_smooth(x) = 0 if x <= -1,
// He_smooth(x) = 1 if x >=  1,
// He_smooth(x) = 0.5 + [15/16].x - [10/16].x^3 + [3/16].x^5 if -1 < x < 1.
// Then, He_smooth(x) is twice differentiable everywhere.
// To have He_smooth(x0,x1,x) for any x0 < x1, one uses the formula He_smooth(x0,x1,x) = He_smooth((2.0*x - x1 - x0)/(x1-x0))
//
// Variables and special values:
// -----------------------------
// x0 : point before which He_smooth is 0.
// x1 : point after which He_smooth is 1.
// x : parameter where one wants He(x0,x1,x).
// 0.9375, -0.625, 0.1875, : 15/16, -10/16, 3/16.

double He_smooth (const double x0,const double x1,const double x)
{
  if (x <= x0)
    return 0.0;
  else if (x >= x1)
    return 1.0;
  else
  {
    const double y = (2.0*x - x1 - x0)/(x1 - x0),He_y = 0.5 + y*(0.937 + y*y*(-0.625 + y*y*0.1875));
    
    return He_y;
  }
}
  





// Class splines approximating a function given by a set of points. 
//------------------------------------------------------------------
// SCALAR_TYPE : double or complex.

template <typename SCALAR_TYPE>
class F_interpolated_class
{
public:

  // Constructor:
  // ------------
  //
  // Variables :
  // -----------
  // N : number of points of the discretized function.
  // r_table,F_table : tables of N points having the abscissas and values of F.
  // dr_0,dr_1 : r_table(1) - r_table(0),dr_2 = r_table(2) - r_table(0).
  // dF0 : approximate first derivative of F at the first point.
  // Nm3 : N-3
  // dr_Nm2,dr_Nm3 : r_table(Nm2) - r_table(Nm1),dr_Nm3 = r_table(Nm3) - r_table(Nm1).
  // dF_Nm1 : approximate first derivative of F at the last point.
  
   F_interpolated_class (const unsigned int N_c,
			 const class array<double> &r_table_c,
			 const class array<SCALAR_TYPE> &F_table) 
     : N (N_c),Nm1 (N_c-1),Nm2 (N_c-2),
       one_over_step (1.0/(r_table_c(1) - r_table_c(0))),
       r0 (r_table_c(0)),
       r_Nm1 (r_table_c(N_c-1))
  {
    r_table_ptr = new class array<double> (N);
    class array<double> &r_table = *r_table_ptr;
    for (unsigned int i = 0 ; i < N ; i++) r_table(i) = r_table_c(i);

    equally_spaced = true;
    for (unsigned int i = 1 ; i < Nm1 ; i++)
      if (abs (r_table(i+1) - r_table(i) - 1.0/one_over_step) > precision)
	equally_spaced = false;
    
    const double dr_1 = r_table(1) - r_table(0),dr_2 = r_table(2) - r_table(0);
    const SCALAR_TYPE dF0 = (dr_2*dr_2*F_table(1)-dr_1*dr_1*F_table(2)-(dr_2*dr_2-dr_1*dr_1)*F_table(0))/(dr_2-dr_1)/dr_1/dr_2;

    const unsigned int Nm3 = N-3;
    const double dr_Nm2 = r_table(Nm2) - r_table(Nm1),dr_Nm3 = r_table(Nm3) - r_table(Nm1);
    const SCALAR_TYPE dF_Nm1 = (dr_Nm3*dr_Nm3*F_table(Nm2)-dr_Nm2*dr_Nm2*F_table(Nm3)-(dr_Nm3*dr_Nm3-dr_Nm2*dr_Nm2)*F_table(Nm1))/(dr_Nm3-dr_Nm2)/dr_Nm2/dr_Nm3;

    splines_coeff = new SCALAR_TYPE * [Nm1];
   for (unsigned int i = 0 ; i < Nm1 ; i++) splines_coeff[i] = new SCALAR_TYPE [4];

    splines (F_table,dF0,dF_Nm1);
  }
  
  ~F_interpolated_class (void) 
  {
    delete r_table_ptr;

    for (unsigned int i = 0 ; i < Nm1 ; i++) delete [] splines_coeff[i];
    delete [] splines_coeff;
    }
  
  SCALAR_TYPE operator() (const double r) const;

  const unsigned int N,Nm1,Nm2;     // number of points of the discretized function, N-1, N-2.
  class array<double> *r_table_ptr; // r_table_ptr : pointer on class array containing the values of r_table.

private:
  const double r0,r_Nm1; // r0 = r_table(0), r_Nm1 = r_table(N-1).

  const double one_over_step; // inverse of the length of the interval r_table(1) - r_table(0). 
                                               // It is the inverse of the table step if one has equally spaced abscissas.
  bool equally_spaced;        // true if one has equally spaced abscissas, false if not.

  SCALAR_TYPE **splines_coeff;  // double pointer on SCALAR_TYPE containing the coefficients of the interpolating polynoms.

  unsigned int i_calc (const double r) const;

  void splines (const class array<SCALAR_TYPE> &F_table,const SCALAR_TYPE dF0,const SCALAR_TYPE d2F_Nm1);
};



// Determination of the cubic splines approximating the discrete function F_table. 
//-----------------------------------------------------------------------------
//
// Variables:
// ----------
// u : class array<SCALAR_TYPE>. Work table of dimension N-1.
// d2F_table : class array<SCALAR_TYPE> having the second-derivative of the spline approximation of F in r_table[i].
// term,u: SCALAR_TYPE and class array<SCALAR_TYPE> used to invert the tridiagonal matrix.
// Dr : r_table(i+1)-r_table(i).
// splines_coeff_i : SCALAR_TYPE yable splines_coeff[i] having the 4 spline coefficients for each interval i.

template <typename SCALAR_TYPE>
void F_interpolated_class<SCALAR_TYPE>::splines (const class array<SCALAR_TYPE> &F_table,
						 const SCALAR_TYPE dF0,const SCALAR_TYPE dF_Nm1)
{
  const class array<double> &r_table = *r_table_ptr;

  class array<SCALAR_TYPE> u(Nm1);
  class array<SCALAR_TYPE> d2F_table(N);
 
  d2F_table(0) = -0.5;
  u(0) = (3.0/(r_table[1]-r_table[0]))*((F_table[1]-F_table[0])/(r_table[1]-r_table[0])-dF0);
   
  for (unsigned int i = 1 ; i < Nm1 ; i++)
  {
    const SCALAR_TYPE term = (r_table(i) - r_table(i-1))/(r_table(i+1) - r_table(i-1));
    const SCALAR_TYPE p = term*d2F_table(i-1) + 2.0;
    d2F_table(i) = (term - 1.0)/p;
    u(i) = (F_table(i+1) - F_table(i))/(r_table(i+1) - r_table(i)) - (F_table(i) - F_table(i-1))/(r_table(i) - r_table(i-1));
    u(i) = (6.0*u(i)/(r_table(i+1) - r_table(i-1)) - term*u(i-1))/p;
  }
 
   const SCALAR_TYPE qn=0.5;
   const SCALAR_TYPE un=(3.0/(r_table[Nm1]-r_table[Nm2]))*(dF_Nm1-(F_table[Nm1]-F_table[Nm2])/(r_table[Nm1]-r_table[Nm2]));
   d2F_table[Nm1]=(un-qn*u[Nm2])/(qn*d2F_table[Nm2]+1.0);
  
  for (unsigned int i = Nm2 ; i <= Nm2 ; i--)
    d2F_table(i) = d2F_table(i)*d2F_table(i+1) + u(i);
  
  for (unsigned int i = 0 ; i < Nm1 ; i++)
  { 
    const double Dr = r_table(i+1) - r_table(i),one_over_Dr = 1.0/Dr;
    SCALAR_TYPE *const splines_coeff_i = splines_coeff[i];

    splines_coeff_i[3] = (d2F_table(i+1) - d2F_table(i))*one_over_Dr*0.166666666666667;
    splines_coeff_i[2] = d2F_table(i)*0.5;
    splines_coeff_i[1] = (F_table(i+1) - F_table(i))*one_over_Dr - (splines_coeff_i[3]*Dr + splines_coeff_i[2])*Dr;
    splines_coeff_i[0] = F_table(i);
  }
}





// Determination of the interval in which the value r is to calculate F(r).
// ------------------------------------------------------------------------
// Knowing r, the interval i in which r is is the largest integer close to (r-r0)*one_over_step if one has equally spaced abscissas,
// except if (r-r0)*one_over_step = N-1 : the interval index is then N-2 because of boundary effects.
// If one does not have equally spaced abscissas, one has to find r with a bisection method.
// 
// Variables :
// -----------
// i_of_r : index of the interval where r is.
// i_deb, i_end : indices with which the bisection method is done.
// x_deb : r-r_table(i_deb).
// i_mid : (i_end + i_deb)/2
// x_mid : r-r_table(i_mud).

template <typename SCALAR_TYPE>
unsigned int F_interpolated_class<SCALAR_TYPE>::i_calc (const double r) const
{
  if (r > r_Nm1) return Nm2;
  if (r < r0) return 0;

  if (equally_spaced == true) 
  {
    const unsigned int i_of_r = make_uns_int (floor ((r-r0)*one_over_step));
    
    if (i_of_r == Nm1) 
      return Nm2; 
    else 
      return i_of_r;
  }

  const class array<double> &r_table = *r_table_ptr;

  unsigned int i_deb = 0,i_end = Nm1;
  double x_deb = r-r0;

  do
  {
    const unsigned int i_mid = (i_end + i_deb)/2;
    const double x_mid = r-r_table(i_mid);
    
    if (SIGN (x_deb) != SIGN (x_mid))
      i_end = i_mid;
    else 
    {
      i_deb = i_mid;
      x_deb = x_mid;
    }
  }
  while (i_end - i_deb > 1);

  return i_deb;
}






// Calculation of the function F in r by spline interpolation.
// -----------------------------------------------------------
// 
// Variables :
// -----------
// i,x : x=r-r_table(i), with i the interval index : r >= r_table(i), and r < r_table(i+1) if r < r_table(N-1). 
// If r < r0 or r > r_Nm1, one uses respectively the spline of i=0 or i=N-2. 
// Fr : F(r) calculated by splines interpolation.

template <typename SCALAR_TYPE>
SCALAR_TYPE F_interpolated_class<SCALAR_TYPE>::operator() (const double r) const
{
  const class array<double> &r_table = *r_table_ptr;
  const unsigned int i = i_calc (r);

  const SCALAR_TYPE Fr = poly_eval<SCALAR_TYPE> (3,splines_coeff[i],r-r_table(i));

  return Fr;  
}
