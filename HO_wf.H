void GL_function(double x1,double x2, class array <double> *  & x, class array <double>*  & w,int n)
{
  
  double Precision=1e-10;
  
  x= new class array<double> (n);
  w= new class array<double> (n);
  
  int m,j,i;     
  double z1,z,xm,xl,pp,p3,p2,p1;
  
  m=(n+1)/2;
  xm=0.5*(x1+x2);
  xl=0.5*(x2-x1);
  
  for (int i=1;i<=m;i++)
    {
      z=cos(M_PI*(i-0.25)/(n+0.5));
      
      do 
        {
          p1=1.0;
          p2=0.0;
          
          for (j=1;j<=n;j++)
            {
              p3=p2;
              p2=p1;
              p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;      
            }
          
          pp=n*(z*p1-p2)/(z*z-1.0);
          z1=z;
          z=z1-p1/pp;
          
        }
      
      while (abs(z-z1)>Precision);
      
      (*x)(i-1)=xm-xl*z;
      (*x)(n-i)=xm+xl*z;
      (*w)(i-1)=2.0*xl/( (1.0-z*z)*pp*pp);
      (*w)(n-i)=(*w)(i-1);      
    }
}



struct HO_wf
{
  int n,l,N_pts;
  double b,R,energy;


  class array<double> *u_wf;
  class array<double> *u_wf_hp;
  class array<double> *u_wf_fixed_step;
  

  double R_wf_in_0;
  double dR_wf_in_0;
  double div_R_wf_in_0;


  double u (const double b,const int n,const int l,const double r);
  double du (const double b,const int n,const int l,const double r);
  //  double d2u (const class HO_class &HO,
  //      const int n,const int l,const double r);
  // double energy (const double hbar_omega,const int n,const int l);

  void constr_u_wf(class array <double> *x);
  
  void constr_u_wf_hp(class array <double> *x);
  
  void constr_u_wf_fixed_step();
 
  void calc_in_0();
  

  double calc_in_r(double r);

};





// Calculation of the Laguerre polynom L(n,l+0.5,x) for HO u(n,l) . 
// --------------------------------------------------------------- 
//
// Variables:
// ----------
// n,l : quantum numbers of the HO state.
// x : variable of the Laguerre polynom.
// y : (-x)^i/i! in the loop.
// binom : binomial coefficient.

double Laguerre (const int n,const int l,const double x)
{
  double L = 0.0, y = 1;
                     
  for (int i = 0 ; i <= n; i++)
  {
    double binom = 1.0;
    for (int j = 1 ; j <= n-i ; j++)
      { 
	binom *= (n + l + 0.5 - j + 1)/j;
      }


    L += binom*y;
 
    y *= (-x)/(i+1.0);
  }

  return L;
}
 
                       

double Laguerre_exp_m_r2_over_2b2 (const int n,const int l,const double x,double b)
{
  double L = 0.0, y = 1;
                     
  for (int i = 0 ; i <= n; i++)
  {
    double binom = 1.0;
    for (int j = 1 ; j <= n-i ; j++)
      binom *= (n + l + 0.5 - j + 1)/j;

    L += binom*y*exp(-x*0.5);
 
    y *= (-x)/(i+1.0);
  }

  return L;
}
 

// Calculation of the Laguerre polynom derivative dL/dx(n,l+0.5,x) for HO du(n,l). 
// -------------------------------------------------------------------------------
// 
// Variables:
// ----------
// n,l : quantum numbers of the HO state.
// x : variable of the Laguerre polynom.
// y : (-x)^(i-1)/(i-1)! in the loop.
// binom : binomial coefficient.

double dLaguerre (const int n,const int l,const double x)
{
  double dL = 0.0, y = 1;
                     
  for (int i = 1 ; i <= n ; i++)
  {
    double binom = 1.0;

    for (int j = 1 ; j <= n-i ; j++)
      binom *= (n + l + 0.5 - j + 1)/j;

    dL -= binom*y;

    y *= (-x)/static_cast<double> (i);
  }

  return dL;
}






// Calculation of u(n,l) (HO).     
// ---------------------------
// 
// Variables:
// ----------
// b : parameter of the HO potential (fm).
// n,l : quantum numbers of the HO state.
// r : radius (fm).

// double HO_wf::u (const double b,const int n,const int l,const double r)     
// {
//   if (r == 0.0) return 0.0;

//   const double norm_HO = sqrt (2.0*exp (lgamma (n+1.0)-lgamma (n+l+1.5))/b);
//   const double rho = r/b, rho_square = rho*rho;

//   return (norm_HO*exp (-0.5*rho_square + (l+1.0)*log (rho))*Laguerre (n,l,rho_square));
// }


void HO_wf::constr_u_wf (class array <double> *x )
{

  u_wf=new class array<double> (N_pts);

  for (int i=0;i<=N_pts-1;i++)
    {
      double r=(*x)(i);

      const double norm_HO = sqrt (2.0*exp (lgamma (n+1.0)-lgamma (n+l+1.5))/b);
      const double rho = r/b, rho_square = rho*rho;

      (*u_wf)(i)=norm_HO*exp (-0.5*rho_square + (l+1.0)*log (rho))*Laguerre (n,l,rho_square);
    }
}



double HO_wf::calc_in_r ( double r )
{
  double u_wf_in_r;

  const double norm_HO = sqrt (2.0*exp (lgamma (n+1.0)-lgamma (n+l+1.5))/b);
  const double rho = r/b, rho_square = rho*rho;
  
  u_wf_in_r=norm_HO*exp (-0.5*rho_square + (l+1.0)*log (rho))*Laguerre (n,l,rho_square);
  
  return  u_wf_in_r;
}



void HO_wf::constr_u_wf_hp (class array <double> *x )
{

  u_wf_hp=new class array<double> (N_pts);

  for (int i=0;i<=N_pts-1;i++)
    {
      double r=(*x)(i);

      const double norm_HO = sqrt (2.0*exp (lgamma (n+1.0)-lgamma (n+l+1.5))/b);
      const double rho = r/b, rho_square = rho*rho;

      (*u_wf_hp)(i)=norm_HO*exp ( (l+1.0)*log (rho))*Laguerre_exp_m_r2_over_2b2(n,l,rho_square,b);
    }

}



void HO_wf::constr_u_wf_fixed_step ()
{
  u_wf_fixed_step=new class array<double> (N_pts);


  double step=R/(N_pts-1);

  for (int i=0;i<=N_pts-1;i++)
    {
      double r=i*step;

      const double norm_HO = sqrt (2.0*exp (lgamma (n+1.0)-lgamma (n+l+1.5))/b);
      const double rho = r/b, rho_square = rho*rho;

      (*u_wf_fixed_step)(i)=norm_HO*exp (-0.5*rho_square + (l+1.0)*log (rho))*Laguerre (n,l,rho_square);
    }
}


// Calculation of du(n,l)/dr (HO). 
// -------------------------------
// 
// Variables:
// ----------
// b : parameter of the HO potential (fm).
// n,l : quantum numbers of the HO state.
// r : radius (fm).

double HO_wf::du (const double b,const int n,const int l,const double r)     
{  
  if ((r == 0) && (l > 0)) return 0.0;

  const double norm_HO = sqrt (2.0*exp (lgamma (n+1.0)-lgamma (n+l+1.5))/b)/b;
  const double rho = r/b, rho_square = rho*rho;

  const double Laguerre_term = Laguerre (n,l,rho_square)*(l+1.0-rho_square) + 2.0*rho_square*dLaguerre (n,l,rho_square);

  if ((r == 0) && (l == 0)) return (norm_HO*Laguerre_term);

  return (norm_HO*exp (-0.5*rho_square + l*log (rho))*Laguerre_term);
}      




void HO_wf::calc_in_0()
{
  R_wf_in_0=du(b,n,l,0.0);

   double r=0.0;

  const double norm_HO = sqrt (2.0*exp (lgamma (n+1.0)-lgamma (n+l+1.5))/b);
  const double rho = r/b, rho_square = rho*rho;
  
  div_R_wf_in_0=norm_HO*exp (-0.5*rho_square)*Laguerre (n,l,rho_square)/pow(b,l+1);
  
  if (l>1)
    {
      div_R_wf_in_0=0.0;
    }


  double div_du_in_0=0.0;

  if (l==1)
    {
      const double norm_HO_bis = sqrt (2.0*exp (lgamma (n+1.0)-lgamma (n+l+1.5))/b)/b;
      
      const double Laguerre_term = Laguerre (n,l,rho_square)*(l+1.0-rho_square) + 2.0*rho_square*dLaguerre (n,l,rho_square);

      div_du_in_0=norm_HO_bis*exp(-0.5*rho_square)*Laguerre_term/pow(b,l);
    }

  
  dR_wf_in_0=div_du_in_0-div_R_wf_in_0;

}





// Calculation of d2u(n,l)/dr^2 (HO). 
// ---------------------------------
// 
// Variables:
// ----------
// HO : HO potential class.
// b : parameter of the HO potential (fm).
// kin_fact : 2.mu/(hbar^2) in the HO equation.
// hbar_omega : hw of the HO potential.
// n,l : quantum numbers of the HO state.
// r : radius (fm).
/*
double HO_wf::d2u (const class HO_class &HO,
		   const int n,const int l,const double r)
{
  const double b = HO.b,kin_fact = HO.kin_fact,hbar_omega = HO.hbar_omega;

  return (l*(l+1)/(r*r) + r*r/(b*b*b*b) - kin_fact*energy (hbar_omega,n,l))*u (b,n,l,r);
}


*/




// Calculation of the energy of u(n,l) (HO). 
// -----------------------------------------
// 
// Variables:
// ----------
// hbar_omega : hw of the HO potential.
// // n,l : quantum numbers of the HO state.

// double HO_wf::energy (const double hbar_omega,const int n,const int l)
// {
//   return (hbar_omega*(2*n + l + 1.5));
// }

