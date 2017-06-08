/**************************************
 * Weichuan Li Non_Local differential equation 
 *************************************/
#include <iostream>
#include <fenv.h>   // enable floating point trap
#include <fstream>
#define MAXN 1000
#define MAXX0 1000

#include <iomanip>
#include <complex>
#include <cmath>
#include <stdio.h>
#include <gsl/gsl_integration.h>
#include<gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf.h>

#include <armadillo>
#include <ctime>// include this header 
#include <gsl/gsl_sf_coulomb.h> 
#include <fenv.h>   // enable floating point trap
#include <fstream>
#define MAXN 1000
#define MAXX0 1000
#include <complex>
#include <cmath>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <math.h>       /* pow */

using std::cerr;
using std::endl;
#include <fstream>
#include <gsl/gsl_sf_coulomb.h> 
using std::ofstream;
#include <cstdlib> // for exit function
//#include <conio.h>
using namespace arma;
using namespace std;

//using namespace Eigen;
using namespace std;
using std::cerr;
using std::endl;

using std::ofstream;
#include <cstdlib> // for exit function
//#include <conio.h>
using namespace arma;
using namespace std;

const double hbarc = 197.3269718;
const double mass_unit=931.494;
const double c_constant=(pow(hbarc,2)/(2.0*mass_unit));

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


   
void gauleg(const double x1, const double x2, std::vector<double> &x, std::vector<double> &w)
{
	const double EPS=1.0e-8;
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
		x[i]=(xm-xl*z+1.0)/2.0;
		x[n-1-i]=(xm+xl*z+1.0)/2.0;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n-1-i]=w[i];
	}
}



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
 
double central_potential(int l,double R,double u)
{


double constant;

constant=c_constant/u;

return constant*((l*(l+1.0))/(R*R));
}



double Vd_potential(double Vd, double R,int n,int z,double rvd,double avd)

{
    double final,rvvd;
    rvvd=rvd*pow(n+z,1.0/3.0);
    return (-4.0*Vd*exp((R-rvvd)/avd))/pow((1.0+exp((R-rvvd)/avd)),2.0) ;
 


}




double V_sp_potential(double Vso,double Rso,double aso,double r,double j,int l,int n,int z)

{

   double Vspin,Rrso;
   Rrso=Rso*pow(n+z,1.0/3.0);
   Vspin=(j*(j+1.0)-l*(l+1.0)-0.5*(0.5+1.0))*(-Vso/(aso*r))*(pow(hbarc/139.6,2))*(exp((r-Rrso)/aso))/(pow((1.0+exp((r-Rrso)/aso)),2.0));

   return Vspin;


}

double W_sp_potential(double Wso,double Rwso,double awso,double r,double j,int l, int n, int z)

{
 
   double Wspin,Rwwso;
   Rwwso=Rwso*pow(n+z,1.0/3.0);
   Wspin=(j*(j+1.0)-l*(l+1.0)-0.5*(0.5+1.0))*(-Wso/(awso*r))*(pow(hbarc/139.6,2))*(exp((r-Rwwso)/awso))/(pow((1.0+exp((r-Rwwso)/awso)),2.0));
 
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

std::complex<double> local_potential(double R,int l,double j,int n,int z1,int z2,double Vv,double rv,double av,double Wv,double rwv,double awv,double Vd,double rvd,double avd,double Wd,double rwd,double awd,double Vso,double Rso,double aso,double Wso,double Rwso,double awso)


{  

    const   complex<double> i(0.0,1.0);    

    return V_potential(Vv,R,n,z1,rv,av)+1.0*i*Wv_potential(Wv,R,n,z1,rwv,awv)+Vd_potential(Vd,R,n,z1,rvd,avd)+1.0*i*Wd_potential(Wd,R,n,z1,rwd,awd)+V_sp_potential(Vso,Rso,aso,R,j,l,n,z1)+1.0*i*W_sp_potential(Wso,Rwso,awso,R,j,l,n,z1)+columb_potential(z1,z2,R,n);

}


std::complex < double >** local_potential_matrix(int N,double a,int l,double jj, int n,int z1,int z2,double Vv,double rv,double av,double Wv,double rwv,double awv,double Vd,double rvd,double avd,double Wd,double rwd,double awd,double Vso,double Rso,double aso,double Wso,double Rwso,double awso)

{   
    std::vector <double> N_basis, N_weight;

    const   complex<double> i(0.0,1.0);    
    double constant;
    double part1,part2,part3,part4;
   
    
    
     gauleg(-1, 1, N_basis,N_weight);
    std::complex < double >** table = new std::complex < double >*[N];
    for(int i = 0; i < N; i++) 

       {
        table[i] = new std::complex < double >[N]; 
    for(int j = 0; j < N; j++)
      {

     if (i==j)
{      
       
       table[i][j] = local_potential(a*N_basis[i],l,jj,n,z1,z2,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,Vso,Rso,aso,Wso,Rwso,awso);

}
     else

{    
   
     
      table[i][j]=0;

}
    }  
    }
    return table;
}

struct my_f_params {int a; double b;}; 
double f1 (double x, void * p) { //Deriviative part of the integral
        struct my_f_params * params = (struct my_f_params *)p;
        int l = (params->a);
        double mu = (params->b);
        
        return exp(x*mu)*Legendre(l,x);
       }
double integration(int l, double mu){
        const   complex<double> i(0.0,1.0); 
        complex <double> c;
        gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

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



std::complex < double > non_local(double r_minus,double r_plus,int l,double jj,int n,int z1,int z2, double beta,double Vv1,double rv1,double av1,double Wv1,double rwv1,double awv1,double Vd1,double rvd1,double avd1,double Wd1,double rwd1,double awd1,double Vso1,double Rso1,double aso1,double Wso1,double Rwso1,double awso1)
    {
    double r_average,zz,c_m;
    const double PI = 3.141592653589793;
    r_average=(r_plus+r_minus)/2.0;
     
    complex < double> c_n,c,total;
    zz=2.0*r_minus*r_plus/(beta*beta);
    const   complex<double> i(0.0,1.0); 
    c_n=2.0*pow(1.0*i,l)*zz;
    c=1.0/(2.0*pow(1.0*i,l));
    complex <double> V_potential_average;
    V_potential_average=local_potential(r_average,l,jj,n,z1,z2,Vv1,rv1,av1,Wv1,rwv1,awv1,Vd1,rvd1,avd1,Wd1,rwd1,awd1,Vso1,Rso1,aso1,Wso1,Rwso1,awso1);
    if ( abs(zz)<=700)
      {  complex <double> part_b,part_a;
     
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




complex<double>** total_potential_matrix(int N,double a,int l,double jj, int n,int z1,int z2,double Vv,double rv,double av,double Wv,double rwv,double awv,double Vd,double rvd,double avd,double Wd,double rwd,double awd,double Vso,double Rso,double aso,double Wso,double Rwso,double awso,double beta,double Vv1,double rv1,double av1,double Wv1,double rwv1,double awv1,double Vd1,double rvd1,double avd1,double Wd1,double rwd1,double awd1,double Vso1,double Rso1,double aso1,double Wso1,double Rwso1,double awso1)

{   
    //std::vector <double> N_basis;
    std::vector <double> N_basis(N);

    std::vector <double> N_weight(N);

     gauleg(-1, 1, N_basis,N_weight);
    
    const   complex<double> i(0.0,1.0);    
    double constant;

  
    complex <double> total,part1,part2;
    std::complex < double >** table = new std::complex < double >*[N];
 for(int i = 0; i < N; i++) 

       {
    table[i] = new std::complex < double >[N];
    
    for(int j = 0; j < N; j++)
      {  

        part2=non_local(a*N_basis[i],a*N_basis[j],l,jj,n,z1,z2,beta,Vv1,rv1,av1,Wv1,rwv1,awv1,Vd1,rvd1,avd1,Wd1,rwd1,awd1,Vso1,Rso1,aso1,Wso1,Rwso1,awso1);
      
        if (i==j)
        {
       part1=local_potential(a*N_basis[i],l,jj,n,z1,z2,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,Vso,Rso,aso,Wso,Rwso,awso);
      
       
       table[i][j] =a*sqrt(N_weight[i]*N_weight[j]/4.0)*part2+part1;
        }
       else
      {
   //    part2=non_local(a*N_basis[i],a*N_basis[j],l,j,n,z1,z2,beta,Vv1,rv1,av1,Wv1,rwv1,awv1,Vd1,rvd1,avd1,Wd1,rwd1,awd1,Vso1,Rso1,aso1,Wso1,Rwso1,awso1);
      
       table[i][j] =a*sqrt(N_weight[i]*N_weight[j]/4.0)*part2;
       }
    
       }
     
       }
 
   
    return table;  

    
}






double Lagrange_function(int i,double r,double a,int N)
      {
     // std::vector <double> N_basis;
      double xuu;
     // N_basis=basis(N,1e-10);
       std::vector <double> N_basis(N);

    std::vector <double> N_weight(N);

     gauleg(-1, 1, N_basis,N_weight);

      if (r==a*N_basis[i-1])
       {
       xuu= pow(-1.0,N+i)*sqrt(a*N_basis[i-1]*(1.0-N_basis[i-1]));
        }
       else 
      {
        xuu= pow(-1.0,N+i)*r/(a*N_basis[i-1])*sqrt(a*N_basis[i-1]*(1.0-N_basis[i-1]))*Legendre(N,2.0*r/a-1.0)/(r-a*N_basis[i-1]);

       }  
     //  cout<<"i"<<" "<<N_basis[0];
       return xuu;
      }



std::complex<double>** Kinetic_Matrix(int N,int l, double u,double a,double B,double E)

{   
 
    const   complex<double> i(0.0,1.0);    
    double constant;
     std::vector <double> N_basis(N);

    std::vector <double> N_weight(N);

     gauleg(-1, 1, N_basis,N_weight);
    double part1,part2,part3,part4;
	cout.precision(15);
  //  cout<<"Constant"<<c_constant<<endl;
    constant=c_constant/u;
      
    std::complex < double >** table = new std::complex < double >*[N];
    for(int i = 0; i < N; i++) 

       {
    table[i] = new std::complex < double >[N];
    for(int j = 0; j < N; j++)
      {
    
     if (i==j)
{      
    
       part1=(4.0*N*N+4.0*N+3.0)*N_basis[i]*(1.0-N_basis[i])-6.0*N_basis[i]+1.0;
       part2=3.0*a*a*(N_basis[i]*N_basis[i])*((1.0-N_basis[i])*(1.0-N_basis[i]));
     
       table[i][j] = constant*(part1/part2)+central_potential(l,a*N_basis[i],u)+(-constant*B/a*Lagrange_function(i+1,a,a,N)*Lagrange_function(j+1,a,a,N))-E;
     // table[i][j] = constant;

}
     else

{    
    
      
      part3=pow(-1.0,i+j)/(a*a*sqrt(N_basis[i]*N_basis[j]*(1.0-N_basis[i])*(1.0-N_basis[j])));
      part4=(N*N*1.0+N*1.0+1.0+(N_basis[i]+N_basis[j]-2.0*N_basis[i]*N_basis[j])/((N_basis[i]-N_basis[j])*(N_basis[i]-N_basis[j]))-1.0/(1.0-N_basis[i])-1.0/(1.0-N_basis[j]));

     table[i][j]=constant*part3*part4+(-constant*B/a*Lagrange_function(i+1,a,a,N)*Lagrange_function(j+1,a,a,N));
    // table[i][j]=constant;

}
    } 
    }

     

  /*  for(int i = 0; i < N; i++) 

       {
    
    for(int j = 0; j < N; j++)
      {
      cout<<table[i][j]<<" ";
      }

     cout<<endl;
     }

*/

   return table; 
}




std::complex<double>** Total_Matrix(int N,int l, double u,double a,double B,double E)

{   
 
    const   complex<double> i(0.0,1.0);    
    double constant;
     std::vector <double> N_basis(N);

    std::vector <double> N_weight(N);

     gauleg(-1, 1, N_basis,N_weight);
    double part1,part2,part3,part4;
	cout.precision(15);
  //  cout<<"Constant"<<c_constant<<endl;
    constant=c_constant/u;
      
    std::complex < double >** table = new std::complex < double >*[N];
    for(int i = 0; i < N; i++) 

       {
    table[i] = new std::complex < double >[N];
    for(int j = 0; j < N; j++)
      {
    
     if (i==j)
{      
    
       part1=(4.0*N*N+4.0*N+3.0)*N_basis[i]*(1.0-N_basis[i])-6.0*N_basis[i]+1.0;
       part2=3.0*a*a*(N_basis[i]*N_basis[i])*((1.0-N_basis[i])*(1.0-N_basis[i]));
     
       table[i][j] = constant*(part1/part2)+central_potential(l,a*N_basis[i],u)+(-constant*B/a*Lagrange_function(i+1,a,a,N)*Lagrange_function(j+1,a,a,N))-E;
     // table[i][j] = constant;

}
     else

{    
    
      
      part3=pow(-1.0,i+j)/(a*a*sqrt(N_basis[i]*N_basis[j]*(1.0-N_basis[i])*(1.0-N_basis[j])));
      part4=(N*N*1.0+N*1.0+1.0+(N_basis[i]+N_basis[j]-2.0*N_basis[i]*N_basis[j])/((N_basis[i]-N_basis[j])*(N_basis[i]-N_basis[j]))-1.0/(1.0-N_basis[i])-1.0/(1.0-N_basis[j]));

     table[i][j]=constant*part3*part4+(-constant*B/a*Lagrange_function(i+1,a,a,N)*Lagrange_function(j+1,a,a,N));
    // table[i][j]=constant;

}
    } 
    }

    return table;  

}






double delc(int n,int l,double eta){
      if ( n == 0)
      {
        return eta/(1.0+l)-atan(eta/(1.0+l));
      }
      else
      {
        return delc(n-1,l,eta)+eta/(1.0+l+n)-atan(eta/(1.0+l+n));
 
      }
              }

double columb_phase(double eta,int l)
   {
    double  euler=0.5772156649;
    double  phase;
    double  C=0.0;
    if (l>0)
      {
      C=C+1.0/l;
      }
 
    else if (l==0)
      {
       C=0.0;
       }
    phase=delc(901,l,eta)+eta*(C-euler);
    
    return phase;
    }

std::complex <double> I_1( std::complex < double > x,int l)
   { 
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     std::complex <double> part;
     part=sin(x-l*PI/2.0);
     return part;
    }

std::complex <double> O_1( std::complex < double > x,int l)
   { 
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     std::complex <double> part;
     part=cos(x-l*PI/2.0);
     return part;
    }


std::complex<double> H_H_plus(double r,int l,int z1,int z2,double u,double E)
    {

     double SchEqConst,columb;
     double eta,k;
     gsl_sf_result F,G,Fp,Gp; 
     double exp_F,exp_G,t;
      
     
     const   complex<double> i(0.0,1.0);
     SchEqConst=u/c_constant;
     k=sqrt(SchEqConst*E);
   //  eta=(z1*z2*1.43997*931.494*u)/(pow(hbarc,2)*k);
     eta=(z1*z2*1.43997*u/(c_constant*k));
     columb=columb_phase(eta,l);
     int iret = gsl_sf_coulomb_wave_FG_e(eta,k*r,l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
   
     return G.val+1.0*i*F.val;

     }

std::complex<double> H_H_plus_prime(double r,int l,int z1,int z2,double u,double E)
    {

     double SchEqConst,columb;
     double eta,k;
     gsl_sf_result F,G,Fp,Gp; 
     double exp_F,exp_G,t; 
     const   complex<double> i(0.0,1.0);
    
  
     SchEqConst=u/c_constant;
     k=sqrt(SchEqConst*E);
   //  eta=(z1*z2*1.43997*931.494*u)/(pow(hbarc,2)*k);
     eta=(z1*z2*1.43997*u/(c_constant*k));
    
     columb=columb_phase(eta,l);

     int iret = gsl_sf_coulomb_wave_FG_e(eta,k*r,l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
   
  
     return k*(Gp.val+1.0*i*Fp.val);

     }





std::complex<double> H_H_minus(double r,int l,int z1,int z2,double u,double E)
    {
    
     double SchEqConst,columb;
     double eta,k;
     gsl_sf_result F,G,Fp,Gp; 
     double exp_F,exp_G,t; 
     const   complex<double> i(0.0,1.0);
   
    
     SchEqConst=u/c_constant;
     k=sqrt(SchEqConst*E);
   //  eta=(z1*z2*1.43997*931.494*u)/(pow(hbarc,2)*k);
     eta=(z1*z2*1.43997*u/(c_constant*k));
    
    
     columb=columb_phase(eta,l);
     int iret = gsl_sf_coulomb_wave_FG_e(eta,k*r,l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);

     return G.val-1.0*i*F.val;

     }

std::complex<double> H_H_minus_prime(double r,int l,int z1,int z2,double u,double E)
    {
   
     double SchEqConst,columb;
     gsl_sf_result F,G,Fp,Gp; 
     double eta,k;
     double exp_F,exp_G,t; 
     const   complex<double> i(0.0,1.0);

  
     SchEqConst=u/c_constant;
     k=sqrt(SchEqConst*E);
   //  eta=(z1*z2*1.43997*931.494*u)/(pow(hbarc,2)*k);
     eta=(z1*z2*1.43997*u/(c_constant*k));
    
    
     columb=columb_phase(eta,l);
     int iret = gsl_sf_coulomb_wave_FG_e(eta,k*r,l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);

     return k*(Gp.val-1.0*i*Fp.val);

     }

  






std::vector<std::complex<double> > Expansion_source_bloch_bound_2(double u, double E,int N,double a,int l,double j,int n,int z1,int z2,double B)
   { 
   //  std::complex <double> u;
     std::complex <double> constant;
     std::complex <double> c;
     double k,eta;
  //   std::cpmplex <double> angle;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793;
     vector <complex<double> > v; 
     
  //   u=m1*m1/(m1+m2);
  //   angle=Phase_shift(m1,m2,E,N,a,l,n,z,Vr,beta,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,B);
     std::complex <double> final;

     
     k=sqrt(-E*u/c_constant);

     
     constant=c_constant/(a*u);
  
  
 
    
     double const h=0.1e-7;
     double df;
    
    
  
     eta=(-z1*z2*1.43997*931.494*u)/(pow(hbarc,2)*k);
     df=(exp(-k*(a+h/2))*pow(2*k*(a+h/2),-eta) - exp(-k*(a-h/2))*pow(2*k*(a-h/2),-eta) )/h;
     c=a*df-B*pow(2*k*a,-eta)*exp((-k)*a) ;
 
       for(int i=0;i<N;i++)
    {  
      
        v.push_back(c*constant*Lagrange_function(i+1,a,a,N));
    }
  
 
 
     return v;
    }



double Hankel_plus(double eta, int l, double k ,double r)
        {
         double final;
         final= gsl_sf_hyperg_U(-eta+l+1,2+2*l,2.0*k*r)*exp(-k*r)*pow(2.0*k*r,1+l);

         return final;

        }


double Hankel_minus(double eta, int l, double k ,double r)
        {
         double final;
         final= gsl_sf_hyperg_U(-eta+l+1,2+2*l,2.0*k*r)*exp(-k*r)*pow(2.0*k*r,1+l);

         return final;

        }


double Whittaker(double eta, int l, double k ,double r)
        {
         double final;
         final= gsl_sf_hyperg_U(-eta+l+1,2+2*l,2.0*k*r)*exp(-k*r)*pow(2.0*k*r,1+l);

         return final;

        }

std::complex <double> Internal_wave_bound(double r, int l, double eta,double k,double a, double a1, std::complex<double> anc1, int N, std::vector<std::complex<double> > &part_a)
      {
    
       
      if (r<=a1) 
      {
      std::complex <double> add_0=0.0; 
      for(int j = 0; j < N; j++)
      {
      add_0=add_0+part_a[j]*Lagrange_function(j+1,r,a,N);
      
      }
      return add_0*anc1;
      }
      else

      {
      return Whittaker(eta,l, k ,r);


      }
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
     cout << radius<<","<<real<<","<<imag<<" ";}
 
     };



void MyClass::setName(double radius1, double real1, double imag1)

{

 vector<MyClass> list;
 radius=radius1;
 real=real1;
 imag=imag1;
 

}

 





std::complex < double >** total_matrix_inverse_2(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j,int n,int z1, int z2,double B)

  { 
   
    std::complex < double >** B_matrix;
    std::complex < double >** E_matrix;
    std::complex < double >** Kinetic_Matrixx;
    std::complex < double >** C_matrixx;
    std::complex <double> part;
    
   Kinetic_Matrixx=Kinetic_Matrix(N,l,u,a,B,E);
 
   std::complex < double >** total_matrix = new std::complex < double >*[N];
   for ( int c = 0 ; c < N ; c++ )
      {
      total_matrix[c] = new std::complex < double >[N]; 
      for ( int d = 0 ; d < N ; d++ )
         {
       
         part=Kinetic_Matrixx[c][d]+matrix[c][d];
         total_matrix[c][d] = part;
      
         }
       }
    cx_mat mat_arma = zeros<cx_mat>(N,N);
    cx_mat inv_mat = zeros<cx_mat>(N,N);
 
    for(int i=N-1;i>=0;i--)
    {
	for(int j=N-1;j>=0;j--)
	{
		mat_arma(i,j) = total_matrix[i][j];
	}
     }
 
    
    inv_mat=mat_arma.i();
   
    for(int i=0;i<N;i++)
    {
	for(int j=0;j<N;j++)
	{     //  cout<<inv_mat(i,j);
		total_matrix[i][j] = inv_mat(i,j);
              //  cout<<total_matrixx[i][j];
   
    }       
    }
   
    complex<double> * mat_ptr = mat_arma.memptr();
   
     return total_matrix;

  }


std::complex <double>  Phase_shift_2(std::complex<double>** matrix,double u,double E,int N,double a,int l, double j,int n,int z1, int z2,double B)
{   
  
    std::complex <double> add;
    std::complex <double> part_b;
    std::complex <double> final;
    vector <complex<double> > part_a;
      std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
    vector <complex<double> > v;
         std::complex <double> S;
     const   complex<double> i(0.0,1.0);
    std::complex < double >** H_matrixx;
 //   H_matrixx=total_matrix_inverse_2(matrix,u, E, N, a, l,j, n, z1,z2,B);
    
    for(int i=1;i<=N;i++)
    {   

        v.push_back(Lagrange_function(i,a,a,N));
    }
 

    for(int i = 0; i < N; i++) 

    {
    std::complex <double> add_0=0.0;  
    for(int j = 0; j < N; j++)
      {
      add_0=add_0+matrix[i][j]*v[j];

      }
    part_a.push_back(add_0);
    }

   
   part_b=0.0;
   for(int r = 0; r < N; r++)
      {
      part_b=part_b+part_a[r]*v[r];

      }
    
   
   R_matrixx=part_b*c_constant/(u*a);
   cout<<"R_matrix"<<" "<<R_matrixx<<endl;
     std::complex <double> part;
     
    

     S=(H_H_minus(a,l,z1,z2,u,E)-a*R_matrixx*H_H_minus_prime(a,l,z1,z2,u,E))/(H_H_plus(a,l,z1,z2,u,E)-a*R_matrixx*H_H_plus_prime(a,l,z1,z2,u,E));
     cout<<"phase shift"<<S<<endl;
    
     phase_shift=log(S)/(1.0*i*2.0);
     cout<<"phase_shift_degree"<<phase_shift<<endl;
     return phase_shift;

  
}




std::complex <double> External_prime_2(double u,double E,int N,double a,int l,double j,int n,int z1 , int z2, double B)
   { 
    // std::complex <double> mu;
     std::complex <double> angle;
     
     double k,eta,t;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 

     
     std::complex <double> final;
     t=c_constant;
  
     k=sqrt(E*u/t);
    
    

    
     eta=(z1*z2*1.43997*u)/(t*k);
  
     final=1.0*i/2.0*(H_H_minus_prime(a,l,z1,z2,u,E)-exp(2.0*i*angle)*H_H_plus_prime(a,l,z1,z2,u,E));
    // final=1.0*i/2.0*(-i*k*exp(-i*(k*a-l*PI/2.0))-exp(2.0*i*angle)*i*k*exp(i*(k*a-l*PI/2.0)));
     return final;
    }

std::complex <double> External_wave_2(std::complex<double> angle,double u,double E,int N,double a,int l, double j, int n,int z1, int z2,double B)
   { 
    // std::complex <double> mu;

     std::complex <double> k;
     std::complex <double> eta;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);
     const double PI = 3.141592653589793; 
     //mu=m1*m1/(m1+m2);
    // angle=Phase_shift_2(matrix, u, E,N,a,l, j, n, z1,z2, B);
     std::complex <double> final;
    
    
     k=sqrt(E*u/c_constant);

     eta=(z1*z2*1.43997*u)/(c_constant*k);
     
   
     final=(cos(angle)+1.0*i*sin(angle))*(cos(angle)*I_1(k*a,l)+sin(angle)*O_1(k*a,l));
       
   
     return final;
    }



std::vector<std::complex<double> > Expansion_source_bloch_2(std::complex<double> angle,double u,double E,int N,double a,int l, double j, int n,int z1, int z2,double B)
   { 
    // std::complex <double> mu;
     std::complex <double> constant;
     std::complex <double> c;
     std::complex <double> k;

     std::complex <double> R_matrixx;

     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793;
     vector <complex<double> > v; 
     
   //  mu=m1*m1/(m1+m2);
   //  angle=Phase_shift_2(matrix,u,E,N,a,l,j,n,z1,z2,B);
     std::complex <double> final;
   
     k=sqrt(E*u/c_constant);
   
     
     constant=c_constant/(a*u);
     c=a*External_prime_2(u, E, N, a, l, j,n, z1,z2, B)-B*1.0*External_wave_2(angle,u, E, N,a,l,j,n,z1,z2, B);
       for(int i=1;i<=N;i++)
    {   

        v.push_back(c*constant*Lagrange_function(i,a,a,N));
    }
  
 
 
     return v;
    }



std::vector<MyClass> Internal_wave_function_2(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j,int n,int z1, int z2, double B, double Nr, double Rmax,std::vector<MyClass>&list)
   { 
    // std::complex <double> mu;
     std::complex <double> constant;
     std::complex <double> c;
     std::complex <double> k;
  
     std::complex <double>** C_matrix;
    
     const   complex<double> i(0.0,1.0);
     vector <complex<double> > b_matrix; 
     vector <complex<double> > a_matrix;
    
     double eta;
     std::complex <double> final,angle;
    
     std::complex < double >** total_matrixxx;

     total_matrixxx= total_matrix_inverse_2(matrix,u, E, N, a, l,j, n, z1,z2, B);
     angle=Phase_shift_2(total_matrixxx, u, E,N,a,l, j, n, z1,z2, B);
     b_matrix=Expansion_source_bloch_2(angle,u, E, N, a, l,j, n, z1,z2, B);

   
    
    
   
   
       
      k=sqrt(E*u/c_constant);
     for(int i = 0; i < N; i++) 

     {
     std::complex <double> add_0=0.0;  
     for(int j = 0; j < N; j++)
      {
      add_0=add_0+total_matrixxx[i][j]*b_matrix[j];

      }
      a_matrix.push_back(add_0);
      }
     
      MyClass object1;
    
      MyClass f12;

     

     
     
       for (double tt = 0.0; tt <= Rmax; tt += Nr)
      {
      
       if (tt<a)
      {
    
    
       std::complex <double> add_0=0.0;  

       for (int j=0; j<N ;j++)
       {
     //  cout<<a_matrix[j]<<" "<<endl;
       add_0=add_0+a_matrix[j]*Lagrange_function(j+1,tt,a,N);
  
       
       }

       f12.setName(tt,add_0.real(),add_0.imag());
      
       list.push_back(f12);
 
       }

       else 
       {
        std::complex <double > rrr,eee,fff;
      
       
      
        rrr=(cos(angle)+1.0*i*sin(angle))*(cos(angle)*I_1(k*tt,l)+sin(angle)*O_1(k*tt,l));
      
        f12.setName(tt,rrr.real(),rrr.imag());
      
        list.push_back(f12);
        }

      
        }
    return list;
}




std::vector<MyClass> Internal_wave_function_bound_2(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j,int n,int z1,int z2, double B, double Nr, double Rmax,std::vector<MyClass>&list)
   { 
    // std::complex <double> u;
     std::complex <double> constant;
     std::complex <double> c;
     double k,eta;
      
     double a1;
     std::complex <double> anc1,anc;
     //std::complex <double> anc1;
     std::complex <double> norm;
     vector <complex<double> > b_matrix; 
     vector <complex<double> > a_matrix;
     const   complex<double> i(0.0,1.0); 
     
     std::complex <double> final;
    
     std::complex < double >** total_matrixxx;


     total_matrixxx= total_matrix_inverse_2(matrix,u, E, N, a, l,j, n, z1,z2, B);

     b_matrix=Expansion_source_bloch_bound_2(u,E,N,a,l,j,n,z1,z2,B);
      
     a1=7;

     
   
    
     for(int ii = 0; ii < N; ii++) 

     {
     std::complex <double> add_0=0.0;  
     for(int jj = 0; jj < N; jj++)
      {
      
      add_0=add_0+total_matrixxx[ii][jj]*b_matrix[jj];
      
      }
      
      a_matrix.push_back(add_0);
      }
    k=sqrt(-E*u/c_constant);
    cout<<"k"<<k<<endl;
    eta=(-z1*z2*1.43997*mass_unit*u)/(pow(hbarc,2)*k);
    cout<<"eta"<<eta<<endl;
    cout<<Internal_wave_bound(2.1,l,eta, k,a,a1,anc1,N,a_matrix)<<endl;
    std::complex <double> sum=0.0;
    for (int iii = 0; iii < N; iii++) 

    {
          //cout<<a_matrix[iii]<<endl;
          sum=sum+a_matrix[iii]*Lagrange_function(iii+1,a1,a,N);

     }
   
    anc1=Whittaker(eta,l,k ,a1)/sum;
    norm=0.0;
    cout<<"anc1"<<anc1;
    cout<<"sum"<<sum<<endl;
    cout<<"norm"<<norm<<endl;
    cout<<Internal_wave_bound(2.1,l,eta, k,a,a1,anc1,N,a_matrix)<<endl;
    for (double iiii = 0.0; iiii <= 28.5; iiii += 0.1)
  
    {
        norm=norm+0.1*std::pow(Internal_wave_bound(iiii,l,eta, k,a,a1,anc1,N,a_matrix),2);
     }
    cout<<"norm"<<norm;
    anc=1.0/sqrt(norm);
    cout<<"anc"<<anc;
    
    
      MyClass object1;
    
      MyClass f12;
      
    
       for (double tt = 0.1; tt <= Rmax; tt += Nr)
      {

       f12.setName(tt,(anc*Internal_wave_bound(tt,l,eta,k,a1,a,anc1, N, a_matrix)).real(),(anc*Internal_wave_bound(tt,l,eta,k,a1,a,anc1, N, a_matrix)).imag());
      // f12.display();    
       list.push_back(f12);

 
       }
  //  list.display();
 
    return list;
   

}







struct Final_return
{
  complex<double>** Green_function;
  std::vector<MyClass> wave_function;
};



Final_return non_local_wavefunction_matrix_in(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j,int n,int z1,int z2, double B, double Nr, double Rmax,Final_return m)
{   

  
   m.Green_function=total_matrix_inverse_2(matrix,u,E,N,a,l,j,n,z1,z2,B);
   if (E>0)
   {
   E=E*u;
   m.wave_function=Internal_wave_function_2(matrix,u,E,N,a,l,j,n,z1,z2, B, Nr, Rmax,m.wave_function);
    }
   
   else

   {

  m.wave_function=Internal_wave_function_bound_2(matrix,u,E,N,a,l,j,n,z1,z2, B, Nr, Rmax,m.wave_function);
   }
   return m;
  
}

Final_return non_local_wavefunction_matrixx_in(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j,int n,int z1,int z2, double B, double Nr, double Rmax,Final_return m)
{


  ofstream outdata; // outdata is like cin
  outdata.open("examplee.txt"); // opens the file
  m.Green_function=non_local_wavefunction_matrix_in(matrix,u,E,N,a,l,j,n,z1,z2, B,Nr,Rmax,m).Green_function;
  m.wave_function=non_local_wavefunction_matrix_in(matrix,u,E,N,a,l,j,n,z1,z2, B,Nr,Rmax,m).wave_function;
/* for( int j=0; j<m.wave_function.size(); j++) {
            std::cout << "WAVE"<<" "<<m.wave_function[j].radius<<" "<<m.wave_function[j].real<<" "<<m.wave_function[j].imag << "  "<<endl;
        }
     

*/
  for (int i=0; i<m.wave_function.size();i++)
 {
 
 outdata<< m.wave_function[i].radius<<"  "<<m.wave_function[i].real<<"  "<< m.wave_function[i].imag<<endl;
 }
 outdata.close();
  return m;
  
}


int main()
{ 
    int N;
    ofstream outtime; // outdata is like cin
    outtime.open("example_time_time_NONLOCAL1.txt"); // opens the file
    for ( N=48; N<50; N=N+2)
     {double T1;
      int start_s=clock(); 
      Final_return b;
      Final_return b2;
      double mu;
      complex <double>** matrix_x;
      mu=40.0/41.0;
      matrix_x=total_potential_matrix(N,15,1,0.5,20,20, 0,50.0,1.25,0.25,15.0,1.25,0.25,0,1.25,0.25,0,1.25,0.25,0.0,1.25,0.25,0.0,1.25,0.25,0.85,50.0,1.25,0.25,15.0,1.22,0.47,0,1.25,0.25,0,1.25,0.25,0.0,1.25,0.25,0,1.25,0.25);
      non_local_wavefunction_matrixx_in(matrix_x,mu,25,N,15,1,0.5,20,20,0,0,0.1,15,b);
         
      int stop_s=clock();
  
      T1=(stop_s-start_s)/double(CLOCKS_PER_SEC);
      cout<< N <<"  " <<T1<<endl;
      outtime<< N<<"  "<<abs(T1)<<"  "<<endl;
  

     }

 

 
 return 0;



}

