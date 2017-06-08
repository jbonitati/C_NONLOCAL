/**************************************
 * Weichuan Li Non_Local differential equation 
 *************************************/
#include <iostream>
#include <fenv.h>   // enable floating point trap
#include <fstream>
#define MAXN 1000
#define MAXX0 1000
#include <complex>
#include <cmath>
#include <vector>
#include <math.h>       /* pow */
//#include <array>
//#include <vectors.H>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
//#include <Eigen/Dense>
//#include <arrays2.H>
//#include <vectors.H>
#include <iostream>
#include <fstream>
#include <assert.h>  
#include <math.h>
#include "matrix.h"
//using namespace Eigen;
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


double Bisect(int n,double a,double b,double xacc)
{
 double mid;
 double fa,fmid;

 while(b-a > xacc)
 {
   mid = a + (b-a)/2.0;
   fa = Legendre(n,a);
   fmid = Legendre(n,mid);
   if(fa*fmid > 0.0)
     a = mid;
   else
     b = mid;
 }

 return mid;
}

std::vector<double> roots(int N, double accuracy)
{

 std::vector <double> v;
 std::vector <double> u;
 
 
 int i,j,r,s;
 double xa,xb; // end point of intervals bracketing the root
 //double * x0 = new double [MAXN][MAXX0];
 double x0[MAXN][MAXX0];
// std::vector<std::double> x;




 for(i=1;i<=N;i++) // polynomial order i=1,...,N
 {

  
  
  for(j=0;j<i;j++) // jth root of order i, j=0,1,...,i-1 
  {
    if(j==0)
     xa = -1.0;
    else
     xa = x0[i-1][j-1];

    if(j==i-1)
     xb = 1.0;
    else
     xb = x0[i-1][j];

    x0[i][j]=Bisect(i,xa,xb,accuracy);
 
    v.push_back(x0[i][j]);
   
  }
  

  
 }
//std::cout << "mylist stores " << v.size() << " numbers.\n";

for(r=v.size()-N;r<=v.size()-1;r++)
    {   

        u.push_back(v[r]);
   
      //  cout<<v[r];
       // cout<<" ";
    }

//for(s=0;s<=u.size()-1;s++)
//    {   

       // u.push_back(v[s]);
   
//        cout<<u[s];
//        cout<<" ";
//    }

return u;
}




std::vector<double> basis(int N, double accuracy)
{

 std::vector <double> v;
 std::vector <double> u;
 
 
 int i,j,r,s;
 double xa,xb; // end point of intervals bracketing the root
 //double * x0 = new double [MAXN][MAXX0];
 double x0[MAXN][MAXX0];
// std::vector<std::double> x;




 for(i=1;i<=N;i++) // polynomial order i=1,...,N
 {

  
  
  for(j=0;j<i;j++) // jth root of order i, j=0,1,...,i-1 
  {
    if(j==0)
     xa = -1.0;
    else
     xa = x0[i-1][j-1];

    if(j==i-1)
     xb = 1.0;
    else
     xb = x0[i-1][j];

    x0[i][j]=Bisect(i,xa,xb,accuracy);
 
    v.push_back(x0[i][j]);
   
  }
  

  
 }

for(r=v.size()-N;r<=v.size()-1;r++)
    {   

        u.push_back((v[r]+1.0)/2.0);
    }


return u;
}




std::complex<double>** S_Matrix(int N,int l, double m1,double m2,double a)

{   
    std::vector <double> N_basis;
    double u;
    const   complex<double> i(0.0,1.0);    
    double constant;
    double part1,part2,part3,part4;
    u=m1*m2/(m1+m2);
    constant=20.736/u;
    N_basis=basis(N,1e-10);
   // std::vector <double> N_basis;
   // std::complex < double >** table = new std::complex < double >*[N];
    std::complex < double >** table = new std::complex < double >*[N];
    for(int i = 0; i < N; i++) 

       {
    table[i] = new std::complex < double >[N];
    for(int j = 0; j < N; j++)
      {

     if (i==j)
{      
       cout<<"basis"<<" ";
       
       cout<<N_basis[i]<<" ";
       cout<<"\n";
       part1=(4.0*N*N+4.0*N+3.0)*N_basis[i]*(1.0-N_basis[i])-6.0*N_basis[i]+1.0;
       part2=3.0*a*a*(N_basis[i]*N_basis[i])*((1.0-N_basis[i])*(1.0-N_basis[i]));
     
       table[i][j] = constant*(part1/part2);

}
     else

{    
    //  part3=sqrt(0.0);
      
      part3=pow(-1.0,i+j)/(a*a*sqrt(N_basis[i]*N_basis[j]*(1.0-N_basis[i])*(1.0-N_basis[j])));
      part4=(N*N*1.0+N*1.0+1.0+(N_basis[i]+N_basis[j]-2.0*N_basis[i]*N_basis[j])/((N_basis[i]-N_basis[j])*(N_basis[i]-N_basis[j]))-1.0/(1.0-N_basis[i])-1.0/(1.0-N_basis[j]));

      table[i][j]=constant*part3*part4;

}
    }// sample set value;    
    }
    return table;
}






double Vd_potential(double Vd, double R,int n,int z,double rvd,double avd)

{
    double final,rvvd;
    rvvd=rvd*pow(n+z,1.0/3.0);
    return (-4.0*Vd*exp((R-rvvd)/avd))/(1.0+exp((R-rvvd)/avd))*(2.0) ;
 


}


double  Wd_potential(double Wd,double R,int n,int z,double rwd,double awd)

{
    double rwwd;
    rwwd=rwd*pow(n+z,1.0/3.0);
 
    return (-4.0*Wd*exp((R-rwd)/awd))/(1.0+exp((R-rwd)/awd))*(2.0) ;


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



double central_potential(int l, double R,double m1,double m2)
{   

    double u, constant;

    u=m1*m2/(m1+m2);
    
  
    constant=20.90/u;

    return constant*((l*(l+1.0))/(R*R));

}


std::complex<double> local_potential(double m1,double m2,double R,int l,int n,int z,double Vv,double rv,double av,double Wv,double rwv,double awv,double Vd,double rvd,double avd,double Wd,double rwd,double awd)


{  

    const   complex<double> i(0.0,1.0);    

    return V_potential(Vv,R,n,z,rv,av)+i*Wv_potential(Wv,R,n,z,rwv,awv)+Vd_potential(Vd,R,n,z,rvd,avd)+i*Wd_potential(Wd,R,n,z,rwd,awd)+central_potential(l,R,m1,m2);

}

std::complex < double >** energy_mmatrix(int N, double E, double m1,double m2)
  


   {   
    std::vector <double> N_basis;
    std::complex < double >** table = new std::complex < double >*[N];
    //double** table = new double*[N];
    for(int i = 0; i < N; i++) 

       {
        table[i] = new std::complex < double >[N];
    for(int j = 0; j < N; j++)
      {

     if (i==j)
{      
       
     
       table[i][j] = -E;

}
     else

{    
    
      

      table[i][j]=0;

}
    }
    }
    return table;
}

//def total_matrix(m1,m2,E,N,a,l,n,z,Vr,beta,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,B):

std::complex < double >** local_potential_matrix(double m1,double m2,int N,double a,int l,int n,int z,double Vv,double rv,double av,double Wv,double rwv,double awv,double Vd,double rvd,double avd,double Wd,double rwd,double awd)

{   
    std::vector <double> N_basis;
    double u;
    const   complex<double> i(0.0,1.0);    
    double constant;
    double part1,part2,part3,part4;
    u=m1*m2/(m1+m2);
    constant=20.736/u;
    N_basis=basis(N,1e-10);
    std::complex < double >** table = new std::complex < double >*[N];
    for(int i = 0; i < N; i++) 

       {
        table[i] = new std::complex < double >[N]; 
    for(int j = 0; j < N; j++)
      {

     if (i==j)
{      
       
       table[i][j] = local_potential(m1,m2,a*N_basis[i],l,n,z,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd);

}
     else

{    
    //  part3=sqrt(0.0);
      
     
      table[i][j]=0;

}
    }// sample set value;    
    }
    return table;
}



double f (int l, void * params) {
  double x = *(double *) params;
  double mu = *(double *) params;
  double f = exp(x*mu)*Legendre(l, mu);
  return f;
}


double ff(int l,double mu,double x)
{


return exp(x*mu)*Legendre(l, mu);


}



double Lagrange_function(int i,double r,double a,int N)
      {
      std::vector <double> N_basis;
      double xuu;
      N_basis=basis(N,1e-10);
      if (r==a*N_basis[i-1])
       {
       xuu= pow(-1.0,N+i)*sqrt(a*N_basis[i]*(1.0-N_basis[i]));
        }
       else 
      {
        xuu= pow(-1.0,N+i)*r/(a*N_basis[i])*sqrt(a*N_basis[i]*(1.0-N_basis[i]))*Legendre(N,2.0*r/a-1.0)/(r-a*N_basis[i]);

       }  
       return xuu;
      }



std::complex < double >** local_B_matrix(double m1,double m2,int N,double a,double B)

{   
    std::vector <double> N_basis;
    double u;
    double constant;
    
    u=m1*m2/(m1+m2);
    constant=20.736/u;
    N_basis=basis(N,1e-10);
    std::complex < double >** table = new std::complex < double >*[N];
    for(int i = 0; i < N; i++) 

       {
        table[i] = new std::complex < double >[N]; 
    for(int j = 0; j < N; j++)
      {

     
       
       table[i][j] = -constant*B/a*Lagrange_function(i+1,a,a,N)*Lagrange_function(j+1,a,a,N);


    }// sample set value;    
    }
    return table;
}




std::complex < double >** total_matrix(double m1,double m2,double E,int N,double a,int l,int n,int z,double Vv,double rv,double av,double Wv,double rwv,double awv,double Vd,double rvd,double avd,double Wd,double rwd,double awd,int B)

  {
    std::complex < double >** total_matrix;
    std::complex < double >** local_matrix;
    std::complex < double >** B_matrix;
    std::complex < double >** E_matrix;
  //  total_matrix= local_B_matrix(m1,m2,N,a,B)+local_potential_matrix(m1,m2,N,a,l,n,z,Vv, rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd);
   local_matrix=local_potential_matrix(m1,m2,N,a,l,n,z,Vv, rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd);
   B_matrix=local_B_matrix(m1,m2,N,a,B);
   E_matrix=energy_mmatrix(N,E,m1,m2);
    
   for ( int c = 0 ; c < N ; c++ )
      {
      for ( int d = 0 ; d < N ; d++ )
         {
         total_matrix[c][d] = local_matrix[c][d] + B_matrix[c][d]+E_matrix[c][d];
         }
       }
    return total_matrix;
  }





std::complex < double >** total_matrix_inverse(double m1,double m2,double E,int N,double a,int l,int n,int z,double Vv,double rv,double av,double Wv,double rwv,double awv,double Vd,double rvd,double avd,double Wd,double rwd,double awd,int B)

  {
    std::complex < double >** total_matrixx;
    std::complex < double >** total_matrixx_inverse;
    
    total_matrixx=total_matrix(m1, m2, E, N, a, l, n, z, Vv, rv, av, Wv, rwv, awv, Vd, rvd, avd, Wd, rwd, awd, B);
   // total_matrixx_inverse=total_matrixx.inverse();
   // total_matrixx_inverse=MatrixInversion(total_matrixx, N, total_matrixx_inverse);
    return inv(total_matrixx);
  }


std::complex <double> H_minus( std::complex < double > k,double r,int l)
   { 
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     std::complex <double> part;
     part=-1.0*i*(k*r-0.5*l*PI);
     return exp(part);
    }



std::complex <double> H_minus_prime( std::complex < double > k,double r,int l)
   { 
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     std::complex <double> part;
     part=-1.0*i*(k*r-0.5*l*PI);
     return -1.0*i*k*exp(part);
    }

std::complex <double> H_plus( std::complex < double > k,double r,int l)
   { 
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     std::complex <double> part;
     part=1.0*i*(k*r-0.5*l*PI);
     return exp(part);
    }



   
std::complex <double> H_plus_prime( std::complex < double > k,double r,int l)
   { 
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     std::complex <double> part;
     part=1.0*i*(k*r-0.5*l*PI);
     return 1.0*i*k*exp(part);
    }




std::complex <double>  R_matrix(double m1,double m2, double E,int N,double a,int l,int n,int z,double Vr, double beta,double Vv,double rv,double av, double Wv,double rwv,double awv,double Vd,double rvd,double avd,double Wd,double rwd,double awd,int B)
{   
    std::complex <double> u;
    std::complex <double> add;
    std::complex <double> part_b;
    std::complex <double> final;
    vector <complex<double> > part_a;
 
    vector <complex<double> > v;
    u=m1*m2/(m1+m2);
    std::complex < double >** H_matrixx;
    H_matrixx=total_matrix_inverse(m1, m2, E, N, a, l, n, z, Vv, rv, av, Wv, rwv, awv, Vd, rvd, avd, Wd, rwd, awd, B);
 
    for(int i=1;i<=N;i++)
    {   

        v.push_back(Lagrange_function(i,a,a,N));
    }
 

    for(int i = 0; i < N; i++) 

    {
    std::complex <double> add_0=0.0;  
    for(int j = 0; j < N; j++)
      {
      add_0=add_0+H_matrixx[i][j]*v[j];

      }
    part_a.push_back(add_0);
    }

   
   part_b=0.0;
   for(int r = 0; r < N; r++)
      {
      part_b=part_b+part_a[r]*v[r];

      }
    
   
   final=part_b*20.903/(u*a);

      return final;

}



std::complex <double> Phase_shift(double m1,double m2, double E,int N,double a,int l,int n,int z,double Vr, double beta,double Vv,double rv,double av, double Wv,double rwv,double awv,double Vd,double rvd,double avd,double Wd,double rwd,double awd,int B)
   { 
     std::complex <double> mu;
     std::complex <double> S;
     std::complex <double> k;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     mu=m1*m1/(m1+m2);
     R_matrixx=R_matrix(m1,m2,E,N,a,l,n,z,Vr,beta,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,B);
     std::complex <double> part;
     if (E>=0)
    {
       k=sqrt(E*mu/20.903);
    }
    else
    {
     
       k=1.0*i*sqrt(-E*mu/20.903);

     }
     S=(H_minus(k,a,l)-a*R_matrixx*H_minus_prime(k,a,l))/(H_plus(k,a,l)-a*R_matrixx*H_plus_prime(k,a,l));
     phase_shift=log(S)/(1.0*i*2.0);
     return phase_shift;
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
   
std::complex <double> I_1_d( std::complex < double > x,int l)
   { 
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     std::complex <double> part;
     part=cos(x-l*PI/2.0);
     return part;
    }


std::complex <double> O_1_d( std::complex < double > x,int l)
   { 
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     std::complex <double> part;
     part=-sin(x-l*PI/2.0);
     return part;
    }

std::complex <double> External_prime(double m1,double m2, double E,int N,double a,int l,int n,int z,double Vr, double beta,double Vv,double rv,double av, double Wv,double rwv,double awv,double Vd,double rvd,double avd,double Wd,double rwd,double awd,int B)
   { 
     std::complex <double> mu;
     std::complex <double> angle;
     std::complex <double> k;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     mu=m1*m1/(m1+m2);
     angle=Phase_shift(m1,m2,E,N,a,l,n,z,Vr,beta,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,B);
     std::complex <double> final;
     if (E>=0)
    {
       k=sqrt(E*mu/20.903);
    }
    else
    {
     
       k=1.0*i*sqrt(-E*mu/20.903);

     }
     
    if (E>=0)
    {
     final=(cos(angle)+1.0*i*sin(angle))*(cos(angle)*k*I_1_d(k*a,l)+sin(angle)*k*O_1_d(k*a,l));
       
    }
    else
    {
     
     final=k*exp(k*a);

     }
     return final;
    }



std::complex <double> External_wave(double m1,double m2, double E,int N,double a,int l,int n,int z,double Vr, double beta,double Vv,double rv,double av, double Wv,double rwv,double awv,double Vd,double rvd,double avd,double Wd,double rwd,double awd,int B)
   { 
     std::complex <double> mu;
     std::complex <double> angle;
     std::complex <double> k;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     mu=m1*m1/(m1+m2);
     angle=Phase_shift(m1,m2,E,N,a,l,n,z,Vr,beta,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,B);
     std::complex <double> final;
     if (E>=0)
    {
       k=sqrt(E*mu/20.903);
    }
    else
    {
     
       k=1.0*i*sqrt(-E*mu/20.903);

     }
     
    if (E>=0)
    {
     final=(cos(angle)+1.0*i*sin(angle))*(cos(angle)*k*I_1(k*a,l)+sin(angle)*k*O_1(k*a,l));
       
    }
    else
    {
     
     final=k*exp(k*a);

     }
     return final;
    }


std::vector<std::complex<double> > Expansion_source_bloch(double m1,double m2, double E,int N,double a,int l,int n,int z,double Vr, double beta,double Vv,double rv,double av, double Wv,double rwv,double awv,double Vd,double rvd,double avd,double Wd,double rwd,double awd,int B)
   { 
     std::complex <double> mu;
     std::complex <double> constant;
     std::complex <double> c;
     std::complex <double> k;
  //   std::cpmplex <double> angle;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793;
     vector <complex<double> > v; 
     mu=m1*m1/(m1+m2);
  //   angle=Phase_shift(m1,m2,E,N,a,l,n,z,Vr,beta,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,B);
     std::complex <double> final;
     if (E>=0)
    {
       k=sqrt(E*mu/20.903);
    }
    else
    {
     
       k=1.0*i*sqrt(-E*mu/20.903);

     }
     constant=20.903/(a*mu);
     c=a*External_prime(m1,m2,E,N,a,l,n,z,Vr,beta,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,B)-B*1.0*External_wave(m1,m2,E,N,a,l,n,z,Vr,beta,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,B);
       for(int i=1;i<=N;i++)
    {   

        v.push_back(c*constant*Lagrange_function(i,a,a,N));
    }
  
 
 
     return v;
    }









int main()
{ 
 std::vector <double> u;
   int s;
   std::complex < double >** matrixx;
  
   u=roots(4.0,1e-10);
   for(s=0;s<=u.size()-1;s++)
    {   

   
        //cout<<u[s];
        cout<<" ";
     }
  matrixx=S_Matrix(3,0,2.0,2.0,2.0);
  for(int x=0;x<3;x++)  // loop 3 times for three lines
    {
        for(int y=0;y<3;y++)  // loop for the three elements on the line
        {
            cout<<matrixx[x][y];  // display the current element out of the array
            cout<< " ";
        }
    cout<<endl;  // when the inner loop is done, go to a new line
    }
 // matrixx=k_matrix(3,3);
 
 // printMatrix(matrixx,3);


}



