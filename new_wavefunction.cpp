
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
#include <valarray>
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

using namespace std;
typedef complex<double> complejo;
typedef vector<double> vector_dbl;
typedef vector<double> vector_cmpx;
using namespace arma;

const double hbarc = 197.3269718;
const double mass_unit=931.494;
const double c_constant=(pow(hbarc,2)/(2.0*mass_unit));

#define MAX_PTS 6000
#define AMU 931.494
#define HC  197.3269718 // MeV*fm
#define E2HC	0.00729927
#define PI 3.1415926
complejo const I(0., 1.);


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
{ //Deriviative part of the integral
        struct my_f_params * params = (struct my_f_params *)p;
        int l = (params->a);
        double mu = (params->b);
        
        return exp(x*mu)*Legendre(l,x);
}
       
double integration(int l, double mu)
{
        const   complex<double> i(0.0,1.0); 
        complex <double> c;
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
        struct my_f_params;               
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

/*Returns the Coulomb potential for particles with proton numbers z1,z2
 *  at distance R*/
double columb_potential(int z1,int z2,double R,int n)
{
    double Rcoul,Vcoul;
    Rcoul=1.25*pow((n+z1),1.0/3.0);
  
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
 * Vv and Wv are real and imaginary volume potentials
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

std::complex<double> local_potential(double R,int l,double j,int n,
	int z1,int z2,double Vv,double rv,double av,double Wv,
	double rwv,double awv,double Vd,double rvd,double avd,double Wd,
	double rwd,double awd,double Vso,double Rso,double aso,double Wso,
	double Rwso,double awso)
{  
    const   complex<double> i(0.0,1.0);    

    return V_potential(Vv,R,n,z1,rv,av)
		+1.0*i*Wv_potential(Wv,R,n,z1,rwv,awv)
		+Vd_potential(Vd,R,n,z1,rvd,avd)
		+1.0*i*Wd_potential(Wd,R,n,z1,rwd,awd)
		+V_sp_potential(Vso,Rso,aso,R,j,l,n,z1)
		+1.0*i*W_sp_potential(Wso,Rwso,awso,R,j,l,n,z1)
		+columb_potential(z1,z2,R,n);
}

std::complex<double> non_local(double r_minus,double r_plus,int l,
	double jj,int n,int z1,int z2, double beta,double Vv1,double rv1,
	double av1,double Wv1,double rwv1,double awv1,double Vd1,double rvd1,
	double avd1,double Wd1,double rwd1,double awd1,double Vso1,
	double Rso1,double aso1,double Wso1,double Rwso1,double awso1)
{
    double r_average,zz;
  
    r_average=(r_plus+r_minus)/2.0;
     
    complex < double> c_n,c,total;
    zz=2.0*r_minus*r_plus/(beta*beta);
    const   complex<double> i(0.0,1.0); 
    c_n=2.0*pow(1.0*i,l)*zz;
    c=1.0/(2.0*pow(1.0*i,l));
    complex <double> V_potential_average;
    V_potential_average=local_potential(r_average,l,jj,n,z1,z2,Vv1,
		rv1,av1,Wv1,rwv1,awv1,Vd1,rvd1,avd1,Wd1,rwd1,awd1,Vso1,Rso1,
		aso1,Wso1,Rwso1,awso1);
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
   
    return V_potential_average*total;
   // return 0;
}


class MyClass 
{
      public:
      double radius;   
      double real;  
      double imag;
   
   
     void setName (double radius1, double real1, double imag1); // member functions set
     void setValue(double radius, double real, double imag) 
      {
            this->radius=radius;
            this->real = real;
            this->imag = imag;
      }
       
     void display()
     {
		cout << radius<<","<<real<<","<<imag<<" ";
     }
 
};

struct Final_return
{
  complex<double>** Green_function;
  std::vector<MyClass> wave_function;
};


extern void gauleg(const double x1, const double x2, 
	std::vector<double> &x, std::vector<double> &w)
{
	const double EPS=1.0e-14;
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	int n=x.size();
	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=0;i<m;i++) {
		z=cos(3.141592654*(i+0.75)/(n+0.5));
		do {
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
		x[i]=xm-xl*z;
		x[n-1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n-1-i]=w[i];
	}
}

struct distorted_wave 
{
	int id;
	int puntos;
	double radio;
	int l;
	double j;
	int nodos;
	double r[MAX_PTS];
	complejo wf[MAX_PTS];
	double energia;
	float spin;
};

struct lagrange
{
  vector_dbl x;   // Lagrange points
  vector_dbl w;   // Lagrange weights
  vector_dbl r;   // radial grid (length pts), from 0 to a
  vector_dbl rr;  // radial grid (length pts), from 0 to Rmax
  int N;          // size of Lagrange basis
  mat basis;     // Lagrange basis functions (pts x N)
  mat basis_r;
  double a;   // size of box
};
  
void LagrangeBasis(lagrange* lag)
{
  const double EPS=1.0e-10;
  int m,j,i;
  double z1,z,pp,p3,p2,p1;
  m=(lag->N+1)/2;
  for (i=0;i<m;i++) {
    z=cos(PI*(i+0.75)/(lag->N+0.5));
    do {
      p1=1.0;
      p2=0.0;
      for (j=0;j<lag->N;j++) {
	p3=p2;
	p2=p1;
	p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
      }
      pp=lag->N*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while (fabs(z-z1) > EPS);
    lag->x[i]=(-z+1.0)/2.0;
  
    lag->x[lag->N-1-i]=(z+1.0)/2.0;
    lag->w[i]=2./((1.0-z*z)*pp*pp);
    lag->w[lag->N-1-i]=lag->w[i];
  }
  

  for(m=0;m<lag->r.size();m++)
  {
      for(i=1;i<=lag->N;i++)
	  { 
	    if (lag->r[m]==lag->a*lag->x[i-1])
	    {
	      lag->basis(m,i-1) = pow(-1.,lag->N+i)
			*sqrt(lag->a*lag->x[i-1]*(1.-lag->x[i-1]));
	    }
	    else 
	    {
	      lag->basis(m,i-1)= pow(-1.0,lag->N+i)
			*lag->r[m]/(lag->a*lag->x[i-1])
			*sqrt(lag->a*lag->x[i-1]*(1.-lag->x[i-1]))
			*gsl_sf_legendre_Pl(lag->N,2.*lag->r[m]/lag->a-1.)
			/(lag->r[m]-lag->a*lag->x[i-1]);
	    }
	  }
  }

 for(int mm=0;mm<lag->rr.size();mm++)
 { 
	//cout<<"rr.size"<<lag->rr.size()<<endl;
    //cout<<"N"<<lag->N<<endl;
    for(int ii=1;ii<=lag->N;ii++)   
	{  
	  // cout<<" "<<mm<<" "<<ii-1<<" ";
      //  lag->basis_r(mm,ii-1)=0;
	  if (lag->rr[mm]==lag->a*lag->x[ii-1])
	    {
	      lag->basis_r(mm,ii-1)= pow(-1.,lag->N+ii)
			*sqrt(lag->a*lag->x[ii-1]*(1.-lag->x[ii-1]));
	    }
	  else 
	    {
	      lag->basis_r(mm,ii-1)= pow(-1.0,lag->N+ii)
			*lag->rr[mm]/(lag->a*lag->x[ii-1])
			*sqrt(lag->a*lag->x[ii-1]*(1.-lag->x[ii-1]))
			*gsl_sf_legendre_Pl(lag->N,2.*lag->rr[mm]/lag->a-1.)
			/(lag->rr[mm]-lag->a*lag->x[ii-1]);
	    }      
	}
 }
}

complejo interpola2D_cmpxVec(complejo** funcion,vector_dbl r1,vector_dbl r2,
		double posicion1,double posicion2)
{
	int indice1,indice2;
	complejo f11, f12, f21,f22;
	double delta_r1,delta_r2;
	//if ((r1.size()<3)||(r2.size()<3)) Error("Number of points has to be greater than 3 in interpola2D_cmpxVec");
	delta_r1=r1[r1.size()-1]-r1[r1.size()-2];
	indice1 = int(ceil(posicion1/delta_r1)) - 1;
	delta_r2=r2[r2.size()-1]-r2[r2.size()-2];
	indice2 = int(ceil(posicion2/delta_r2)) - 1;
	if (indice1 >r1.size() - 2)
		indice1=r1.size()-2;
	if (indice1 <= 0)
		indice1=1;
	if (indice2 >r2.size() - 2)
		indice2=r2.size()-2;
	if (indice2 <= 0)
		indice2=1;
	f11=funcion[indice1][indice2];
	f12=funcion[indice1][indice2+1];
	f21=funcion[indice1+1][indice2];
	f22=funcion[indice1+1][indice2+1];
	return (f11*(r1[indice1+1]-posicion1)*(r2[indice2+1]-posicion2)
			+f21*(posicion1-r1[indice1])*(r2[indice2+1]-posicion2)
			+f12*(r1[indice1+1]-posicion1)*(posicion2-r2[indice2])
			+f22*(posicion1-r1[indice1])*(posicion2-r2[indice2]))
			/(delta_r1*delta_r2);
}


complejo NLwavefunction(distorted_wave* dw,complejo** v,vector_dbl r1,
	vector_dbl r2, double q1q2, double masa,double radio_max,
	int puntos,double radio_match,ofstream* fp,lagrange* lag)
{       
  double delta_r,hbarx,part3,part4,ri,rj,q,etac,exp_F,exp_G;
  //removed "central" and "delta_a" (unused)
  complejo pot,Rmatrix,Hp,Hm,Hmp,Hpp,S,phase_shift,factor;
  int i,j;
 
  gsl_sf_result F,G,Fp,Gp;
  hbarx=HC*HC/(2.*AMU*masa);
  q=sqrt(dw->energia/hbarx);
  etac=q1q2*masa*E2HC*AMU/(HC*q);
  cx_mat TLmatrix=zeros<cx_mat>(lag->N,lag->N);
  cx_mat Vmatrix=zeros<cx_mat>(lag->N,lag->N);
  cx_mat Gmatrix=zeros<cx_mat>(lag->N,lag->N);
  vec basis_a(lag->N);
  cx_vec c;
  delta_r=radio_max/double(puntos);
  dw->puntos=puntos;
  dw->radio=radio_max;
  basis_a=lag->basis.row(lag->basis.n_rows-1).t(); // vector of length N with the value of each Lagrange function at a
  // Initialization of Coulomb functions at a
   
  gsl_sf_coulomb_wave_FG_e(etac,q*lag->a,dw->l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
  Hp=(G.val+I*F.val);
  Hm=(G.val-I*F.val);
  Hpp=q*(Gp.val+I*Fp.val);
  Hmp=q*(Gp.val-I*Fp.val);

   // Kinetic energy, Bloch term and energy eigenvalue
  for(i=0;i<lag->N;i++)
    {
      TLmatrix(i,i)=hbarx
		*((4.*lag->N*lag->N+4.*lag->N+3.0)*lag->x[i]*(1.-lag->x[i])-6.*lag->x[i]+1.)
		/(3.*lag->a*lag->a*(lag->x[i]*lag->x[i])*((1.-lag->x[i])*(1.-lag->x[i])))-dw->energia;
      for(j=i+1;j<lag->N;j++)
	  {
		  part3=pow(-1.,i+j)/(lag->a*lag->a*sqrt(lag->x[i]*lag->x[j]*(1.-lag->x[i])*(1.-lag->x[j])));
		  part4=(lag->N*lag->N*1.+lag->N*1.0+1.0+(lag->x[i]+lag->x[j]-2.*lag->x[i]*lag->x[j])/
			 ((lag->x[i]-lag->x[j])*(lag->x[i]-lag->x[j]))-1./(1.-lag->x[i])-1./(1.-lag->x[j]));
		  TLmatrix(i,j)=hbarx*part3*part4;
		  TLmatrix(j,i)=TLmatrix(i,j);
	  }
    }
  // Potential (multiplied by ri*rj), including central potential and energy
  for(i=0;i<lag->N;i++)
    { 
      ri=lag->a*lag->x[i];
      Vmatrix(i,i)=hbarx*dw->l*(dw->l+1.)/(ri*ri);  // central potential and energy
      for(j=0;j<lag->N;j++)
	  {          
		  rj=lag->a*lag->x[j];
		  pot=interpola2D_cmpxVec(v,r1,r2,ri,rj);
		  Vmatrix(i,j)+=ri*rj*lag->a*sqrt(lag->w[i]*lag->w[j]/4.)*pot;
	  }
    }
  Gmatrix=inv(TLmatrix+Vmatrix);
 
  Rmatrix=hbarx*dot(basis_a.t(),Gmatrix*basis_a)/lag->a;
  cout<<"part_b"<<endl;
  cout<<"part_b"<<" "<<dot(basis_a.t(),Gmatrix*basis_a)<<endl;
  cout<<"R-matrix: "<<Rmatrix<<endl;
  cout<<"part_b"<<" "<<dot(basis_a.t(),Gmatrix*basis_a)<<endl;
  S=(Hm-lag->a*Rmatrix*Hmp)/(Hp-lag->a*Rmatrix*Hpp);
  cout<<"S_matrix"<<" "<<S<<" "<<endl;
  phase_shift=-I*log(S)/2.;
  cout<<"Phase shift: "<<phase_shift<<endl;
  c=I*hbarx*(Hmp-S*Hpp)*Gmatrix*basis_a/2.;  // c coefficients (see eq. (3.13))
  //abs(c).print(misc4);
  for(i=0;i<lag->basis.n_rows;i++)
    {
      dw->wf[i]=0.;
      ri=(i+1.)*delta_r;
      dw->r[i]=ri;
      for(j=0;j<lag->N;j++)
		{
		  dw->wf[i]+=c(j)*lag->basis(i,j);
		}
    }
   for(i=lag->basis.n_rows;i<dw->puntos;i++)
    {
      dw->wf[i]=0.;
      ri=(i+1.)*delta_r;
      dw->r[i]=ri;
      gsl_sf_coulomb_wave_FG_e(etac,q*ri,dw->l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
      dw->wf[i]=exp(I*(phase_shift))*(cos(phase_shift)*F.val+sin(phase_shift)*G.val);
    }
   gsl_sf_coulomb_wave_FG_e(etac,q*lag->a,dw->l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
  factor=exp(I*(phase_shift))*(cos(phase_shift)*F.val+sin(phase_shift)*G.val)/dw->wf[lag->basis.n_rows];
  for(i=0;i<dw->puntos;i++)
    {
      dw->wf[i]=factor*dw->wf[i];
      *fp<<dw->r[i]<<"   "<<real(dw->wf[i])<<"  "<<imag(dw->wf[i])<<endl;
    }
  //int stop_s=clock();
  //cout<<"Time in NLwavefunction: "<<(stop_s-start_s)/double(CLOCKS_PER_SEC)<<" s"<<endl;
  return phase_shift;
}

complejo NLwavefunction_p(distorted_wave* dw, double q1q2, double masa,double radio_max,
			int puntos,double radio_match,ofstream* fp,lagrange* lag,
			double jj, int n,int z1,int z2,double Vv,double rv,double av,
			double Wv,double rwv,double awv,double Vd,double rvd,double avd,
			double Wd,double rwd,double awd,double Vso,double Rso,
			double aso,double Wso,double Rwso,double awso,double beta,
			double Vv1,double rv1,double av1,double Wv1,double rwv1,
			double awv1,double Vd1,double rvd1,double avd1,double Wd1,
			double rwd1,double awd1,double Vso1,double Rso1,double aso1,
			double Wso1,double Rwso1,double awso1)
{       
  double delta_r,hbarx,part3,part4,ri,rj,q,etac,exp_F,exp_G; 
  //removed "central" and "delta_a" (unused)
  std::complex<double>local,nonlocal;
  complejo pot,Rmatrix,Hp,Hm,Hmp,Hpp,S,phase_shift,factor;
  int i,j;
  int start_s=clock();
  gsl_sf_result F,G,Fp,Gp;
  hbarx=HC*HC/(2.*AMU*masa);
  q=sqrt(dw->energia/hbarx);
  etac=q1q2*masa*E2HC*AMU/(HC*q);
  cx_mat TLmatrix=zeros<cx_mat>(lag->N,lag->N);
  cx_mat Vmatrix=zeros<cx_mat>(lag->N,lag->N);
  cx_mat Gmatrix=zeros<cx_mat>(lag->N,lag->N);
  vec basis_a(lag->N);
  cx_vec c;
  delta_r=radio_max/double(puntos);
  dw->puntos=puntos;
  dw->radio=radio_max;
  basis_a=lag->basis.row(lag->basis.n_rows-1).t(); // vector of length N with the value of each Lagrange function at a
  // Initialization of Coulomb functions at a
 
  gsl_sf_coulomb_wave_FG_e(etac,q*lag->a,dw->l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
  Hp=(G.val+I*F.val);
  Hm=(G.val-I*F.val);
  Hpp=q*(Gp.val+I*Fp.val);
  Hmp=q*(Gp.val-I*Fp.val);

   // Kinetic energy, Bloch term and energy eigenvalue
  for(i=0;i<lag->N;i++)
    { 
      // TLmatrix(i,i)=hbarx*((4.0*lag->N*lag->N+4.0*lag->N+3.0)*lag->x[i]*(1.0-lag->x[i])-6.0*lag->x[i]+1.0)/(3.0*lag->a*lag->a*(lag->x[i]*lag->x[i])*((1.0-lag->x[i])*(1.0-lag->x[i])))-dw->energia;
      TLmatrix(i,i)=hbarx*((4.*lag->N*lag->N+4.*lag->N+3.0)*lag->x[i]*(1.-lag->x[i])-6.*lag->x[i]+1.)/
	(3.*lag->a*lag->a*(lag->x[i]*lag->x[i])*((1.-lag->x[i])*(1.-lag->x[i])))-dw->energia;
      for(j=i+1;j<lag->N;j++)
	{
	  part3=pow(-1.,i+j)/(lag->a*lag->a*sqrt(lag->x[i]*lag->x[j]*(1.-lag->x[i])*(1.-lag->x[j])));
	  part4=(lag->N*lag->N*1.+lag->N*1.0+1.0+(lag->x[i]+lag->x[j]-2.*lag->x[i]*lag->x[j])/
		 ((lag->x[i]-lag->x[j])*(lag->x[i]-lag->x[j]))-1./(1.-lag->x[i])-1./(1.-lag->x[j]));
	  TLmatrix(i,j)=hbarx*part3*part4;
	  TLmatrix(j,i)=TLmatrix(i,j);
	}
    } 
   
   /*for(i=0;i<lag->N;i++)
    {
      
      for(j=0;j<lag->N;j++)
       {
       cout<<" "<<TLmatrix(i,j)<<"  ";
        }

    cout<<endl;
       }
  
*/
  // Potential (multiplied by ri*rj), including central potential and energy
  for(i=0;i<lag->N;i++)
    { 
      ri=lag->a*lag->x[i];
      Vmatrix(i,i)=hbarx*dw->l*(dw->l+1.)/(ri*ri);  // central potential and energy
     
      for(j=0;j<lag->N;j++)
	{          
	  rj=lag->a*lag->x[j];
          nonlocal=non_local(ri,rj,dw->l,jj,n,z1,z2,beta,Vv1,rv1,av1,
			Wv1,rwv1,awv1,Vd1,rvd1,avd1,Wd1,rwd1,awd1,
			Vso1,Rso1,aso1,Wso1,Rwso1,awso1);
	  //pot=interpola2D_cmpxVec(v,r1,r2,ri,rj);
           if (i==j)
          {      
        
        local=local_potential(ri,dw->l,jj,n,z1,z2,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,Vso,Rso,aso,Wso,Rwso,awso);
      //  cout<<" "<<" l"<<dw->l<<" "<<"jj"<<jj<<" "<<"n"<<" "<<n<<" "<<"z1"<<" "<<z1<<" "<<"z2"<<" "<<z2<<" ";
        Vmatrix(i,j)+=lag->a*sqrt(lag->w[i]*lag->w[i]/4.0)*nonlocal+local;
       // cout<<lag->a<<endl;
      //  cout<<" "<<lag->a*sqrt(lag->w[i]*lag->w[i]/4.0)*nonlocal+local<<" ";
          }
           else
          {    
    
        Vmatrix(i,j)+=lag->a*sqrt(lag->w[i]*lag->w[j]/4.0)*nonlocal;
      //  cout<<" "<<lag->a*sqrt(lag->w[i]*lag->w[j]/4.0)*nonlocal<<" ";
           }
	  
	}
    // cout<<endl;
    }
  
 /* for(i=0;i<lag->N;i++)
    {
      
      for(j=0;j<lag->N;j++)
       {
       cout<<TLmatrix+Vmatrix<<"  ";
        }

    cout<<endl;
       }
 */
  Gmatrix=inv(TLmatrix+Vmatrix);
  /* for(i=0;i<lag->N;i++)
    {
      
      for(j=0;j<lag->N;j++)
       {
       cout<<" "<<Gmatrix(i,j)<<"  ";
        }

    cout<<endl;
       }
 
*/
  for(j=0;j<lag->N;j++)
  {
       cout<<" "<<basis_a[j]<<"  ";
  }


  Rmatrix=hbarx*dot(basis_a.t(),Gmatrix*basis_a)/lag->a;
   cout<<"part_b"<<" "<<dot(basis_a.t(),Gmatrix*basis_a)<<endl;
  cout<<"R-matrix: "<<Rmatrix<<endl;
  S=(Hm-lag->a*Rmatrix*Hmp)/(Hp-lag->a*Rmatrix*Hpp);
  
  cout<<"S_matrix"<<" "<<S<<" "<<endl;
  phase_shift=-I*log(S)/2.;
  cout<<"Phase shift: "<<phase_shift<<endl;
  c=I*hbarx*(Hmp-S*Hpp)*Gmatrix*basis_a/2.;  // c coefficients (see eq. (3.13))
  //abs(c).print(misc4);
  for(i=0;i<lag->basis_r.n_rows;i++)
    {
      dw->wf[i]=0.;
      ri=(i+1.)*delta_r;
     // ri=(i+1)*lag->a/lag->N;
      dw->r[i]=ri;
      for(j=0;j<lag->N;j++)
  	{
  	  dw->wf[i]+=c(j)*lag->basis_r(i,j);
          
  	}
     *fp<<dw->r[i]<<"   "<<real(dw->wf[i])<<"  "<<imag(dw->wf[i])<<endl;
    }
  for(i=lag->basis_r.n_rows;i<dw->puntos;i++)
    {
      dw->wf[i]=0.;
      ri=(i+1.)*delta_r;
      dw->r[i]=ri;
      gsl_sf_coulomb_wave_FG_e(etac,q*ri,dw->l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
      dw->wf[i]=exp(I*(phase_shift))*(cos(phase_shift)*F.val+sin(phase_shift)*G.val);
      *fp<<dw->r[i]<<"   "<<real(dw->wf[i])<<"  "<<imag(dw->wf[i])<<endl;
    }
 //  gsl_sf_coulomb_wave_FG_e(etac,q*lag->a,dw->l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
 // factor=exp(I*(phase_shift))*(cos(phase_shift)*F.val+sin(phase_shift)*G.val)/dw->wf[lag->basis.n_rows];
 /* for(i=0;i<dw->puntos;i++)
    {
      //dw->wf[i]=factor*dw->wf[i];
      *fp<<dw->r[i]<<"   "<<real(dw->wf[i])<<"  "<<imag(dw->wf[i])<<endl;
    }
 */
}

int main()
{   
	int start_s=clock();
    ofstream outtime; // outdata is like cin
    outtime.open("wave_function.txt"); // opens the file

      Final_return b;
      Final_return b2;
      //double mu;
     
   boost::property_tree::ptree pt;
   boost::property_tree::ini_parser::read_ini("config.txt", pt);
 int NN;double m1;double m2;double E;double a_size;int l;double jj;
 int n;int z1;int z2;int B; double Nr; double R_max;double Vv;
 double rv;double av;double Wv;double rwv;double awv;double Vd;
 double rvd;double avd;double Wd;double rwd;double awd;double Vso;
 double Rso;double aso;double Wso;double Rwso;double awso;double beta;
 double Vv1;double rv1;double av1;double Wv1;double rwv1;double awv1;
 double Vd1;double rvd1;double avd1;double Wd1;double rwd1;double awd1;
 double Vso1;double Rso1;double aso1;double Wso1;double Rwso1;double awso1;
   
	NN=pt.get<int>("Numerical.NN") ;

	m1=pt.get<double>("Numerical.m1") ;
	m2=pt.get<double>("Numerical.m2") ;

	E=pt.get<double>("Numerical.E") ;
	a_size=pt.get<double>("Numerical.a_size") ;
	l=pt.get<int>("Numerical.l") ;
	jj=pt.get<double>("Numerical.jj");
	n=pt.get<int>("Numerical.n");
	z1=pt.get<int>("Numerical.z1");
	z2=pt.get<int>("Numerical.z2");
	B=pt.get<int>("Numerical.B") ;
	Nr=pt.get<double>("Numerical.Nr");
	R_max=pt.get<double>("Numerical.R_max");


	Vv=pt.get<double>("local.Vv") ;
	rv=pt.get<double>("local.rv") ;
	av=pt.get<double>("local.av") ;
	Wv=pt.get<double>("local.Wv") ;
	rwv=pt.get<double>("local.rwv") ;
	awv=pt.get<double>("local.awv");

	Vd=pt.get<double>("local.Vd");
	rvd=pt.get<double>("local.rvd");
	avd=pt.get<double>("local.avd");

	Wd=pt.get<double>("local.Wd") ;
	rwd=pt.get<double>("local.rwd");
	awd=pt.get<double>("local.awd");

	Vso=pt.get<double>("local.Vso");
	Rso=pt.get<double>("local.Rso");
	aso=pt.get<double>("local.aso");
	Wso=pt.get<double>("local.Wso") ;
	Rwso=pt.get<double>("local.Rwso");
	awso=pt.get<double>("local.awso");

	Vv1=pt.get<double>("Non_local.Vv1") ;
	rv1=pt.get<double>("Non_local.rv1") ;
	av1=pt.get<double>("Non_local.av1") ;
	Wv1=pt.get<double>("Non_local.Wv1") ;
	rwv1=pt.get<double>("Non_local.rwv1") ;
	awv1=pt.get<double>("Non_local.awv1");
	Vd1=pt.get<double>("Non_local.Vd1");
	rvd1=pt.get<double>("Non_local.rvd1");
	avd1=pt.get<double>("Non_local.avd1");
	Wd1=pt.get<double>("Non_local.Wd1") ;
	rwd1=pt.get<double>("Non_local.rwd1");
	awd1=pt.get<double>("Non_local.awd1");
	Vso1=pt.get<double>("Non_local.Vso1");
	Rso1=pt.get<double>("Non_local.Rso1");
	aso1=pt.get<double>("Non_local.aso1");
	Wso1=pt.get<double>("Non_local.Wso1") ;
	Rwso1=pt.get<double>("Non_local.Rwso1");
	awso1=pt.get<double>("Non_local.awso1");

	beta=pt.get<double>("Non_local.beta");
	
	distorted_wave wave;
	lagrange* lag=new lagrange[1]; 

	lag->N=NN;
	lag->a=a_size;
	double q1q2, masa, radio_max, radio_match;
	int puntos;
	q1q2=z1*z2;
	masa=m2/(m1+m2);
	cout<<"masa"<<masa<<endl;
	radio_max=R_max;
	radio_match=a_size;
	puntos=R_max/Nr;

	double step,rn,spin;
	step=double(lag->a/lag->N);
	cout<<"N"<<" "<<lag->N<<" ";
	rn=step;
	
	while(rn<=lag->a)
	{
		  lag->r.push_back(rn);
		  rn+=step;
	}

	/*
	double r=0;
	for(int n=1; n<=lag->N; n++)
	  { 
		r=n*lag->a/lag->N;
		lag->r.push_back(r);
		
	   }
	*/
	cout<<"Size: "<<lag->r.size()<<" Last: "<<lag->r[lag->r.size()-1]<<endl;
	double r_0=0;
	for(int n1=1; n1<=lag->a/Nr; n1++)
	  { 
		r_0=n1*Nr;
		lag->rr.push_back(r_0);
		
	   }
	cout<<"Size: "<<lag->rr.size()<<" Last: "<<lag->rr[lag->rr.size()-1]<<endl;

	/*for(double r_0=0.0; r_0<lag->a; r_0=r_0+Nr)
	  { 
	  
		lag->rr.push_back(r_0);
		
	   }
	cout<<"Size: "<<lag->rr.size()<<" Last: "<<lag->rr[lag->rr.size()-1]<<endl;
	*/
	for(int n1=0;n1<lag->N;n1++)
	{
	  lag->x.push_back(0.);
	  lag->w.push_back(0.);
	}
	lag->basis.zeros(lag->r.size(),lag->N);  
	lag->basis_r.zeros(lag->rr.size(),lag->N); 
	LagrangeBasis(lag);
	distorted_wave* chi=new distorted_wave;

	spin=0.5;
	ofstream fp("neutron_distorted_wave.txt");
	chi->energia=E*masa;
	chi->l=l;
	chi->spin=spin;
	chi->j=jj;
	chi->puntos=puntos;
	NLwavefunction_p(chi,q1q2, masa,radio_max,puntos,radio_match,
		&fp,lag,jj,n,z1,z2,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,
		Vso,Rso,aso,Wso,Rwso,awso,beta,Vv1,rv1, av1,Wv1,rwv1,awv1,Vd1,
		rvd1,avd1, Wd1,rwd1,awd1,Vso1,Rso1,aso1,Wso1,Rwso1,awso1);
	 //non_local_wavefunction_matrixx_in(NN,mu,E,a_size,l,jj,n,z1,z2,B,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,Vso,Rso,aso,Wso,Rwso,awso,beta,Vv1,rv1, av1,Wv1,rwv1,awv1,Vd1,rvd1, avd1, Wd1,rwd1,awd1,Vso1,Rso1,aso1,Wso1,Rwso1,awso1,Nr,R_max,b);
	 int stop_s=clock();
	  cout<<"Time in NLwavefunction: "<<(stop_s-start_s)/double(CLOCKS_PER_SEC)<<" s"<<endl;
	 
	return 0;
}
