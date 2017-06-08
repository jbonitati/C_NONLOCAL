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
#include"lapacke.h"
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


static void __attribute__((constructor)) trapfpe () {
 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}  // program will stop if an invalid number is generated









double Legendre(int n, double t) // return P_{n}(t)
{
 int k;
 double Pk_1,Pk_2,Pk;

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
 double xa,xb; 

 double x0[MAXN][MAXX0];





 for(i=1;i<=N;i++) 
 {

  
  
  for(j=0;j<i;j++)
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


for(r=v.size()-N;r<=v.size();r++)
    {   

        u.push_back(v[r]);
   
   
    }



return u;
}


double central_potential(int l,double R,double u)
{


double constant;

constant=(pow(197.32705,2)/(2.0*931.494))/u;

return constant*((l*(l+1))/(R*R));
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



complex<double> interpola2D_cmpx(complex <double>** funcion,vector<double> r1,vector<double> r2,
double posicion1,double posicion2,int puntos1,int puntos2)
{
int indice1,indice2;
complex<double> f11, f12, f21,f22;
double delta_r1,delta_r2;
//if ((puntos1<3)||(puntos2<3)) Error("The number of points must be greater than 3");
delta_r1=r1[puntos1-1]-r1[puntos1-2];
indice1 = int(ceil(posicion1/delta_r1)) - 1;
delta_r2=r2[puntos2-1]-r2[puntos2-2];
indice2 = int(ceil(posicion2/delta_r2)) - 1;
if (indice1 > puntos1 - 2)
indice1=puntos1-2;
if (indice1 <= 0)
indice1=1;
if (indice2 > puntos2 - 2)
indice2=puntos2-2;
if (indice2 <= 0)
indice2=1;
f11=funcion[indice1][indice2];
f12=funcion[indice1][indice2+1];
f21=funcion[indice1+1][indice2];
f22=funcion[indice1+1][indice2+1];
return (f11*(r1[indice1+1]-posicion1)*(r2[indice2+1]-posicion2)+f21*(posicion1-r1[indice1])*(r2[indice2+1]-posicion2)
+f12*(r1[indice1+1]-posicion1)*(posicion2-r2[indice2])+f22*(posicion1-r1[indice1])
*(posicion2-r2[indice2]))/(delta_r1*delta_r2);
}

vector<double> point_list( )

{
 
 const   complex<double> i(0.0,1.0); 


 int number_of_array;
 int number_of_lines = 0;
 std::string line;
 std::ifstream myfile("pot_p3.txt");

 while (std::getline(myfile, line))
	{
        ++number_of_lines;
	}

 std::cout << "Number of lines in text file: " << number_of_lines<<endl;
 number_of_array=sqrt(number_of_lines);
 cout<<"number of array"<<number_of_array<<endl;
 
 double exam1[number_of_lines];
 double exam2[number_of_lines];
 double exam3[number_of_lines]; 
 double exam4[number_of_lines];

 ifstream infile;   
   std::vector<double> point_array;
 //double* point_array;
 int num = 0; 
 infile.open("pot_p3.txt");
     if(infile.fail()) 
    { 
      cout << "error" << endl; 

    } 
       while(!infile.eof()) 
     {
         infile >> exam1[num]; 
         infile >> exam2[num];
         infile >> exam3[num]; 
          infile >> exam4[num];
         ++num;

        
         // infile >> exam1[num] >> exam2[num] >> exam3[num]; ++num;
      } 
  infile.close(); 
 
  for (int x=0;x<number_of_array;x++)
  {
    point_array.push_back(exam2[x]);
  //  point_array=7.0;

   }



  /* for (int x=0;x<number_of_array;x++)
  {

    cout<<point_array[x]<<" ";
  }
  */
   return point_array;
  
   }
complex<double>** original_matrix()


{ ofstream outdata; // outdata is like cin
  outdata.open("example.txt"); // opens the file
 
      std::complex < double >** matrixx;
   

      const   complex<double> i(0.0,1.0); 
 
      int number_of_array;
      int number_of_lines = 0;
      std::string line;
      std::ifstream myfile("pot_p3.txt");

      while (std::getline(myfile, line))
	{
        ++number_of_lines;
	}

     std::cout << "Number of lines in text file: " << number_of_lines<<endl;
     number_of_array=sqrt(number_of_lines);
     cout<<"number of array"<<number_of_array<<endl;
     std::complex < double >** matrix = new std::complex < double >*[number_of_array];





     double exam1[number_of_lines];
     double exam2[number_of_lines]; 
     double exam3[number_of_lines];
     double exam4[number_of_lines];

     ifstream infile;   

     int num = 0; 
     infile.open("pot_p3.txt");
     if(infile.fail()) 
     { 
      cout << "error" << endl; 
      
     } 
       while(!infile.eof())
      { 
         infile >> exam1[num];
         infile >> exam2[num]; 
         infile >> exam3[num]; 
          infile >> exam4[num];
         ++num; 

      
         // infile >> exam1[num] >> exam2[num] >> exam3[num]; ++num;
      } 
    infile.close(); 
    double (*matrix_real)[number_of_array]= reinterpret_cast<double (*)[number_of_array]>(exam3);
    double (*matrix_imag)[number_of_array]= reinterpret_cast<double (*)[number_of_array]>(exam4);
    for ( int c = 0 ; c <  number_of_array ; c++ )
      {
      matrix[c] = new std::complex < double >[number_of_array]; 
      for ( int d = 0 ; d <  number_of_array ; d++ )
         {
         matrix[c][d]=matrix_real[c][d]+1.0*i*matrix_imag[c][d];

      
         }
       }
  
    return matrix;
 


}



complex<double>** potential_matrix(complex <double>** funcion, vector<double> r1, vector<double> r2,
int N,double a)

{   
    std::vector <double> N_basis;
    ofstream outdata; // outdata is like cin
  outdata.open("example1.txt"); // opens the file
    ofstream outdata2; // outdata is like cin
  outdata2.open("example2.txt"); // opens the file
    const   complex<double> i(0.0,1.0);    
    double constant;
    double part1,part2,part3,part4;
  
    std::vector<std::complex<double> > old;
   // std::vector<std::complex<double> > new;
    N_basis=basis(N,1e-10);
    


   
  for (int i=0; i<150;i++)
 {
 
 outdata<< r1[i]<<"  "<<funcion[i][i].real()<<" "<<funcion[i][i].imag()<<"  "<<endl;
 }
 outdata.close();
  

    std::complex < double >** table = new std::complex < double >*[N];

    for(int i = 0; i < N; i++) 

       {
    table[i] = new std::complex < double >[N];
    
    for(int j = 0; j < N; j++)
      {
    
   
      
       
       table[i][j] =a*N_basis[i]*a*N_basis[j]*float(a)*sqrt(N_weight[i]*N_weight[j]/4.0)*interpola2D_cmpx(funcion,r1,r2,a*N_basis[i],a*N_basis[j],150,150);
       cout<<table[i][j]<<" ";
   

       }
       cout<<endl;
       }

    for (int i=0; i<N;i++)
 {
 
 outdata2<< a*N_basis[i]<<"  "<<table[i][i].real()<<" "<<table[i][i].imag()<<"  "<<endl;
 }
 outdata2.close();
    return table;  

}

complex<double>** central_matrix(int l,double u,int N,double a)

{   
    std::vector <double> N_basis;

    N_basis=basis(N,1e-10);

    std::complex < double >** table = new std::complex < double >*[N];
    for(int i = 0; i < N; i++) 

       {
    table[i] = new std::complex < double >[N];
    for(int j = 0; j < N; j++)
      {
      if (i==j)
       {
   
      
     
       table[i][j] =central_potential(l,a*N_basis[i],u);
    
        }

     else
       {
       table[i][j]=0;
       }
       }
     }
    return table;  

}


complex<double>** optical_matrix( int N,double a)

{   
 vector <double> radial_point;



radial_point=point_list(); 
complex<double>** matrix=original_matrix();


return potential_matrix(matrix,radial_point,radial_point,N,a);


}



std::complex<double>** S_Matrix(int N,int l, double u,double a)

{   
    std::vector <double> N_basis;
   // double u;
    const   complex<double> i(0.0,1.0);    
    double constant;
    double part1,part2,part3,part4;
  //  u=m1*m2/(m1+m2);
    constant=(pow(197.3269718,2)/(2.0*931.494))/u;
    // std::vector <double> N_basis;
    N_basis=basis(N,1e-10);
   // std::complex < double >** table = new std::complex < double >*[N];
    std::complex < double >** table = new std::complex < double >*[N];
    for(int i = 0; i < N; i++) 

       {
    table[i] = new std::complex < double >[N];
    for(int j = 0; j < N; j++)
      {
    
     if (i==j)
{      
    //   cout<<"basis"<<" ";
       
      // cout<<N_basis[i]<<" ";
      // cout<<"\n";
       part1=(4.0*N*N+4.0*N+3.0)*N_basis[i]*(1.0-N_basis[i])-6.0*N_basis[i]+1.0;
       part2=3.0*a*a*(N_basis[i]*N_basis[i])*((1.0-N_basis[i])*(1.0-N_basis[i]));
     
       table[i][j] = constant*(part1/part2);
     // table[i][j] = constant;

}
     else

{    
    
      
      part3=pow(-1.0,i+j)/(a*a*sqrt(N_basis[i]*N_basis[j]*(1.0-N_basis[i])*(1.0-N_basis[j])));
      part4=(N*N*1.0+N*1.0+1.0+(N_basis[i]+N_basis[j]-2.0*N_basis[i]*N_basis[j])/((N_basis[i]-N_basis[j])*(N_basis[i]-N_basis[j]))-1.0/(1.0-N_basis[i])-1.0/(1.0-N_basis[j]));

     table[i][j]=constant*part3*part4;
    // table[i][j]=constant;

}
    } 
    }

    return table;  

}






std::complex < double >** energy_mmatrix(int N, double E, double u)
  


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



double Lagrange_function(int i,double r,double a,int N)
      {
      std::vector <double> N_basis;
      double xuu;
      N_basis=basis(N,1e-10);

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

std::complex < double >** local_B_matrix(double u,int N,double a,double B)

{   
    std::vector <double> N_basis;
   // double u;
    double constant;
    
   // u=m1*m2/(m1+m2);
    constant=(pow(197.3269718,2)/(2.0*931.494))/u;
    N_basis=basis(N,1e-10);
    
    std::complex < double >** table = new std::complex < double >*[N];
    for(int i = 0; i < N; i++) 

       {
        table[i] = new std::complex < double >[N]; 
    for(int j = 0; j < N; j++)
      {
      
        table[i][j] = -constant*B/a*Lagrange_function(i+1,a,a,N)*Lagrange_function(j+1,a,a,N);
      // table[i][j] = -constant*B/a;


}
     
       
      
     //   cout<<-constant*B/a;
       // table[i][j]==-constant*B/a;
    }// sample set value;    
    
    return table;
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

std::complex <double> O_1_2( std::complex < double > x,int l,int z1,int z2,double u,double E)
   { 
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     const double  hbarc=197.3269718;
     double SchEqConst,k,eta,columb;
     std::complex <double> part;
     SchEqConst=(2.0*931.494*u)/(pow(hbarc,2.0));
     k=sqrt(SchEqConst*E);
     eta=(z1*z2*1.43997*931.494*u)/(pow(hbarc,2.0)*k);
     columb=columb_phase(eta,l);
     part=cos(x-l*PI/2.0-eta*log(2.0*x)+columb);
     return part;
    }


std::complex <double> I_1_2( std::complex < double > x,int l,int z1,int z2,double u,double E)
   { 
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     const double  hbarc=197.3269718;
     double SchEqConst,k,eta,columb;
     std::complex <double> part;
     SchEqConst=(2.0*931.494*u)/(pow(hbarc,2.0));
     k=sqrt(SchEqConst*E);
     eta=(z1*z2*1.43997*931.494*u)/(pow(hbarc,2.0)*k);
     columb=columb_phase(eta,l);
     part=sin(x-l*PI/2.0-eta*log(2.0*x)+columb);
     return part;
    }
std::complex<double> H_H_plus(double r,int l,int z1,int z2,double u,double E)
    {
     const double hbarc=197.3269718;
     double SchEqConst,columb;
     double eta,k;
     gsl_sf_result F,G,Fp,Gp; 
     double exp_F,exp_G; 
     const   complex<double> i(0.0,1.0);
     SchEqConst=(2.0*931.494*u)/(pow(hbarc,2.0));
     k=sqrt(SchEqConst*E);
     eta=(z1*z2*1.43997*931.494*u)/(pow(hbarc,2)*k);
     columb=columb_phase(eta,l);
     int iret = gsl_sf_coulomb_wave_FG_e(eta,k*r,l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
   
     return G.val+1.0*i*F.val;

     }

std::complex<double> H_H_plus_prime(double r,int l,int z1,int z2,double u,double E)
    {
     const double hbarc=197.3269718;
     double SchEqConst,columb;
     double eta,k;
     gsl_sf_result F,G,Fp,Gp; 
     double exp_F,exp_G; 
     const   complex<double> i(0.0,1.0);
     SchEqConst=(2.0*931.494*u)/(pow(hbarc,2.0));
     k=sqrt(SchEqConst*E);
     eta=(z1*z2*1.43997*931.494*u)/(pow(hbarc,2)*k);
    // cout<<"eta"<<eta;
     columb=columb_phase(eta,l);
    // cout<<"columb"<<columb<<endl;
     int iret = gsl_sf_coulomb_wave_FG_e(eta,k*r,l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
   
    // double const h=0.1e-7;
    // std::complex <double> df;
    
    // df=( H_H_plus(r+h/2.0,l,z1,z2,u,E) - H_H_plus(r-h/2.0,l,z1,z2,u,E) )/(h);
    // return df;
     return k*(Gp.val+1.0*i*Fp.val);

     }





std::complex<double> H_H_minus(double r,int l,int z1,int z2,double u,double E)
    {
     const double hbarc=197.3269718;
     double SchEqConst,columb;
     double eta,k;
     gsl_sf_result F,G,Fp,Gp; 
     double exp_F,exp_G; 
     const   complex<double> i(0.0,1.0);
     SchEqConst=(2.0*931.494*u)/(pow(hbarc,2.0));
     k=sqrt(SchEqConst*E);
     eta=(z1*z2*1.43997*931.494*u)/(pow(hbarc,2)*k);
     columb=columb_phase(eta,l);
     int iret = gsl_sf_coulomb_wave_FG_e(eta,k*r,l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
    // cout<<"columb"<<columb<<endl;
    // cout<<"g"<<G.val<<endl;
    // cout<<"f"<<F.val<<endl;
     return G.val-1.0*i*F.val;

     }

std::complex<double> H_H_minus_prime(double r,int l,int z1,int z2,double u,double E)
    {
     const double hbarc=197.3269718;
     double SchEqConst,columb;
     gsl_sf_result F,G,Fp,Gp; 
     double eta,k;
     double exp_F,exp_G; 
     const   complex<double> i(0.0,1.0);
     SchEqConst=(2.0*931.494*u)/(pow(hbarc,2.0));
     k=sqrt(SchEqConst*E);
     eta=(z1*z2*1.43997*931.494*u)/(pow(hbarc,2)*k);
     columb=columb_phase(eta,l);
     int iret = gsl_sf_coulomb_wave_FG_e(eta,k*r,l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
   //  double const h=0.1e-7;
   //  std::complex <double> df;
    
    // df=( H_H_minus(r+h/2.0,l,z1,z2,u,E) - H_H_minus(r-h/2.0,l,z1,z2,u,E) )/(h);
    // return df;
     return k*(Gp.val-1.0*i*Fp.val);

     }

  
std::complex <double> I_1_d( std::complex < double > x,int l)
   { 
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     std::complex <double> part;
     part=cos(x-l*PI/2.0);
     return part;
    }

std::complex <double> I_1_d_2( std::complex < double > x,int l,int z1,int z2,double u,double E)
   { 
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     const double  hbarc=197.3269718;
     double SchEqConst,k,eta,columb;
     std::complex <double> part;
     SchEqConst=(2.0*931.494*u)/(pow(hbarc,2.0));
     k=sqrt(SchEqConst*E);
     eta=(z1*z2*1.43997*931.494*u)/(pow(hbarc,2.0)*k);
     columb=columb_phase(eta,l);
     part=cos(x-l*PI/2.0-eta*log(2.0*x)+columb);
     return part;
    }



std::complex <double> O_1_d_2( std::complex < double > x,int l,int z1,int z2,double u,double E)
   { 
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     const double  hbarc=197.3269718;
     double SchEqConst,k,eta,columb;
     std::complex <double> part;
     SchEqConst=(2.0*931.494*u)/(pow(hbarc,2.0));
     k=sqrt(SchEqConst*E);
     eta=(z1*z2*1.43997*931.494*u)/(pow(hbarc,2.0)*k);
     columb=columb_phase(eta,l);
     part=-sin(x-l*PI/2.0-eta*log(2.0*x)+columb);
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

     
     k=sqrt(-E*u/(pow(197.3269718,2)/(2.0*931.494)));

     
     constant=(pow(197.3269718,2)/(2.0*931.494))/(a*u);
     const double  hbarc=197.3269718;
    // std::complex <double> final;
 
    
     double const h=0.1e-7;
     double df;
    
    
  
     eta=(-z1*z2*1.43997*931.494*u)/(pow(hbarc,2)*k);
     df=(exp(-k*(a+h/2))*pow(2*k*(a+h/2),-eta) - exp(-k*(a-h/2))*pow(2*k*(a-h/2),-eta) )/h;
     c=a*df-B*pow(2*k*a,-eta)*exp((-k)*a) ;
  //   c=a*(-k)*exp((-k)*a)-B*exp((-k)*a) ;
       for(int i=0;i<N;i++)
    {  
       // cout<<c*constant*Lagrange_function(i+1,a,a,N);
        v.push_back(c*constant*Lagrange_function(i+1,a,a,N));
    }
  
 
 
     return v;
    }


std::complex <double> Internal_wave(double r, double a, int N, std::vector<std::complex<double> > &part_a)
      {
    
      std::complex <double> add_0=0.0;  

    //  cout<<"part_a";
      for(int j = 0; j < N; j++)
      {
      add_0=add_0+part_a[j]*Lagrange_function(j+1,r,a,N);
    //  cout<<part_a[j];
      }
     

      return add_0;
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

 





std::complex < double >** total_matrix_2(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j,int n,int z1, int z2,double B)

  { 
   // std::complex < double >** total_matrixx;
   // std::complex < double >** local_matrix;
    std::complex < double >** B_matrix;
    std::complex < double >** E_matrix;
    std::complex < double >** S_matrixx;
    std::complex < double >** C_matrixx;
    std::complex <double> part;
  //  total_matrix= local_B_matrix(m1,m2,N,a,B)+local_potential_matrix(m1,m2,N,a,l,n,z,Vv, rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd);
   //local_matrix=local_potential_matrix(m1,m2,N,a,l,n,z,Vv, rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd);
   B_matrix=local_B_matrix(u,N,a,B);
   E_matrix=energy_mmatrix(N,E,u);
   S_matrixx=S_Matrix(N,l,u,a);
   C_matrixx=central_matrix(l,u,N,a);
   std::complex < double >** table = new std::complex < double >*[N];
   for ( int c = 0 ; c < N ; c++ )
      {
      table[c] = new std::complex < double >[N]; 
      for ( int d = 0 ; d < N ; d++ )
         {
         part= matrix[c][d]+ B_matrix[c][d]+E_matrix[c][d]+S_matrixx[c][d]+C_matrixx[c][d];
         table[c][d] = part;
        // cout<<part;
         }
       }
    return table;
  }

std::complex < double >** total_matrix_inverse_2(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j,int n,int z1, int z2,double B)

  {
   
    std::complex < double >** total_matrixx;
    Mat<cx_double> total_matrixx_inverse;
    Mat<cx_double> A;
    total_matrixx=total_matrix_2(matrix,u, E, N, a, l, j, n, z1, z2,B);
   
    cx_mat mat_arma = zeros<cx_mat>(N,N);
    cx_mat inv_mat = zeros<cx_mat>(N,N);
 
    for(int i=N-1;i>=0;i--)
    {
	for(int j=N-1;j>=0;j--)
	{
		mat_arma(i,j) = total_matrixx[i][j];
	}
     }
 
    
    inv_mat=mat_arma.i();
   
    for(int i=0;i<N;i++)
    {
	for(int j=0;j<N;j++)
	{     //  cout<<inv_mat(i,j);
		total_matrixx[i][j] = inv_mat(i,j);
              //  cout<<total_matrixx[i][j];
   
    }       
    }
   
    complex<double> * mat_ptr = mat_arma.memptr();
   
     return total_matrixx;
    }


std::complex <double>  R_matrix_2(std::complex<double>** matrix,double u,double E,int N,double a,int l, double j,int n,int z1, int z2,double B)
{   
   // std::complex <double> u;
    std::complex <double> add;
    std::complex <double> part_b;
    std::complex <double> final;
    vector <complex<double> > part_a;
 
    vector <complex<double> > v;
   // u=m1*m2/(m1+m2);
   
    std::complex < double >** H_matrixx;
    H_matrixx=total_matrix_inverse_2(matrix,u, E, N, a, l,j, n, z1,z2,B);
 
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
    
   
   final=part_b*(pow(197.3269718,2)/(2.0*931.494))/(u*a);
     cout<<"Rmatrix"<<" "<<1.0/final<<endl;

      return final;
   

}

std::complex <double> Phase_shift_2(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j,int n,int z1,int z2,double B)
   { 
     //std::complex <double> mu;
     std::complex <double> S;
     double k;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
    // mu=m1*m1/(m1+m2);
     R_matrixx=R_matrix_2(matrix,u, E, N, a, l, j,n, z1,z2,B);
     std::complex <double> part;
     
    
     k=sqrt(E*u/(pow(197.3269718,2)/(2.0*931.494)));
     S=(H_H_minus(a,l,z1,z2,u,E)-a*R_matrixx*H_H_minus_prime(a,l,z1,z2,u,E))/(H_H_plus(a,l,z1,z2,u,E)-a*R_matrixx*H_H_plus_prime(a,l,z1,z2,u,E));
    
    // S=(H_minus(k,a,l)-a*R_matrixx*H_minus_prime(k,a,l))/(H_plus(k,a,l)-a*R_matrixx*H_plus_prime(k,a,l));
     phase_shift=log(S)/(1.0*i*2.0);
     return phase_shift;
    }


std::complex <double> External_prime_2(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j,int n,int z1 , int z2, double B)
   { 
    // std::complex <double> mu;
     std::complex <double> angle;
     
     double k,eta;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     const double  hbarc=197.3269718;
    // mu=m1*m1/(m1+m2);
     angle=Phase_shift_2(matrix,u, E, N, a, l,j, n, z1,z2,B);
     std::complex <double> final;
     k=sqrt(E*u/(pow(197.3269718,2)/(2.0*931.494)));
    
    
    eta=(z1*z2*1.43997*931.494*u)/(pow(hbarc,2.0)*k);
  
    final=1.0*i/2.0*(H_H_minus_prime(a,l,z1,z2,u,E)-exp(2.0*i*angle)*H_H_plus_prime(a,l,z1,z2,u,E));
    
     return final;
    }

std::complex <double> External_wave_2(std::complex<double>** matrix,double u,double E,int N,double a,int l, double j, int n,int z1, int z2,double B)
   { 
    // std::complex <double> mu;
     std::complex <double> angle;
     std::complex <double> k;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793; 
     //mu=m1*m1/(m1+m2);
     angle=Phase_shift_2(matrix, u, E,N,a,l, j, n, z1,z2, B);
     std::complex <double> final;
    
    
     k=sqrt(E*u/(pow(197.3269718,2)/(2.0*931.494)));
    
     
   
     final=(cos(angle)+1.0*i*sin(angle))*(cos(angle)*I_1(k*a,l)+sin(angle)*O_1(k*a,l));
       
   
     return final;
    }



std::complex <double> Internal_wave_funtion_1(double Nr,double Rmax)
   { 
 
  
     // MyClass object1;
      vector<MyClass> list;
      MyClass t_123;
   
      cout<<"Rmax"<<Rmax<<" "<<"Nr"<<" "<<Nr<<" ";
      for (double tt = 0.1; tt <= Rmax; tt += Nr)
      {
       //cout<<t;
       t_123.setName(tt,tt,tt);
       t_123.display();    
 
 
       }
 
  
  }
std::vector<std::complex<double> > Expansion_source_bloch_2(std::complex<double>** matrix,double u,double E,int N,double a,int l, double j, int n,int z1, int z2,double B)
   { 
    // std::complex <double> mu;
     std::complex <double> constant;
     std::complex <double> c;
     std::complex <double> k;
  //   std::cpmplex <double> angle;
     std::complex <double> R_matrixx;
     std::complex <double> phase_shift;
     const   complex<double> i(0.0,1.0);  
     const double PI = 3.141592653589793;
     vector <complex<double> > v; 
     
   //  mu=m1*m1/(m1+m2);
  //   angle=Phase_shift(m1,m2,E,N,a,l,n,z,Vr,beta,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,B);
     std::complex <double> final;
   
     k=sqrt(E*u/(pow(197.3269718,2)/(2.0*931.494)));
   
     
     constant=(pow(197.3269718,2)/(2.0*931.494))/(a*u);
     c=a*External_prime_2(matrix,u, E, N, a, l, j,n, z1,z2, B)-B*1.0*External_wave_2(matrix,u, E, N,a,l,j,n,z1,z2, B);
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
    
     
     vector <complex<double> > b_matrix; 
     vector <complex<double> > a_matrix;
    
  
     std::complex <double> final;
    
     std::complex < double >** total_matrixxx;
    // vector<MyClass> list;
     total_matrixxx= total_matrix_inverse_2(matrix,u, E, N, a, l,j, n, z1,z2, B);
     //  total_matrix_inverse(m1,m2,E,N,a,l,n,z,Vr,beta,Vv,rv,av,Wv,rwv,awv,Vd,rvd,avd,Wd,rwd,awd,B);
     b_matrix=Expansion_source_bloch_2(matrix,u, E, N, a, l,j, n, z1,z2, B);

    
     for(int i = 0; i < N; i++) 

     {
     std::complex <double> add_0=0.0;  
     for(int j = 0; j < N; j++)
      {
      add_0=add_0+total_matrixxx[i][j]*b_matrix[j];

      }
      a_matrix.push_back(add_0);
      }
     // cout<<"inter"<< Internal_wave(0.1,a, N, a_matrix).real();
      MyClass object1;
     // vector<MyClass> list;
      MyClass f12;
   
     // cout<<"Rmax"<<Rmax<<" "<<"Nr"<<" "<<Nr<<" ";
       for (double tt = 0.1; tt <= Rmax; tt += Nr)
      {
      // cout<<"t";
       f12.setName(tt,Internal_wave(tt,a, N, a_matrix).real(),Internal_wave(tt,a, N, a_matrix).imag());
      // f12.display();    
       list.push_back(f12);

 
       }
  //  list.display();
 
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
     const double  hbarc=197.3269718;
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
    k=sqrt(-E*u/(pow(197.3269718,2)/(2.0*931.494)));
    cout<<"k"<<k<<endl;
    eta=(-z1*z2*1.43997*931.494*u)/(pow(hbarc,2)*k);
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
      
     // cout<<"Rmax"<<Rmax<<" "<<"Nr"<<" "<<Nr<<" ";
       for (double tt = 0.1; tt <= Rmax; tt += Nr)
      {
     //  cout<<t;
       f12.setName(tt,(anc*Internal_wave_bound(tt,l,eta,k,a1,a,anc1, N, a_matrix)).real(),(anc*Internal_wave_bound(tt,l,eta,k,a1,a,anc1, N, a_matrix)).imag());
      // f12.display();    
       list.push_back(f12);

 
       }
  //  list.display();
 
    return list;
   

}




std::vector<MyClass> Internal_wave_function_general_2(std::complex<double>** matrix,double u,double E,int N,double a,int l,double j,int n,int z1,int z2, double B, double Nr, double Rmax,std::vector<MyClass>&list)
   { 
    
     if (E>=0)
    {  E=E*u;
       list=Internal_wave_function_2(matrix,u,E,N,a,l,j,n,z1,z2, B, Nr, Rmax,list);
    }
    else
    {
       
       list=Internal_wave_function_bound_2(matrix,u,E,N,a,l,j,n,z1,z2, B, Nr, Rmax,list);

     }
    
    return list;
    }


std::vector<MyClass> Internal_wave_function_general_3(std::complex<double>** function, vector<double> r1,vector <double> r2,double u,double E,int N,double a,int l,double j,int n,int z1,int z2, double B, double Nr, double Rmax,std::vector<MyClass>&list)
   { 
     
      complex<double>**matrix;
 //  complex<double>**Green_function;
   matrix=potential_matrix(function,r1, r2,N,a);

    
     if (E>=0)
    {  E=E*u;
       list=Internal_wave_function_2(matrix,u,E,N,a,l,j,n,z1,z2, B, Nr, Rmax,list);
    }
    else
    {
       
       list=Internal_wave_function_bound_2(matrix,u,E,N,a,l,j,n,z1,z2, B, Nr, Rmax,list);

     }
    
    return list;
    }

struct Final_return
{
  complex<double>** Green_function;
  std::vector<MyClass> wave_function;
};


Final_return non_local_wavefunction_matrix(std::complex<double>** function,vector<double> r1,vector<double> r2,double u,double E,int N,double a,int l,double j,int n,int z1,int z2, double B, double Nr, double Rmax,Final_return m)
{   

  // std::vector<MyClass> wave_function;
   complex<double>**matrix;
 //  complex<double>**Green_function;
   matrix=potential_matrix(function,r1, r2,N,a);

   m.Green_function=total_matrix_inverse_2(matrix,u,E,N,a,l,j,n,z1,z2,B);
   
   m.wave_function=Internal_wave_function_general_3(function,r1,r2,u,E,N,a,l,j,n,z1,z2,B,Nr,Rmax,m.wave_function);
   return m;
  
}





Final_return non_local_wavefunction_matrixx(double u,double E,int N,double a,int l,double j,int n,int z1,int z2, double B, double Nr, double Rmax,Final_return m)
{

  std::complex<double>**original;
  std::vector<double> list;
  
  original=original_matrix();
  list=point_list();
 
  ofstream outdata; // outdata is like cin
  outdata.open("example.txt"); // opens the file
  m.Green_function=non_local_wavefunction_matrix(original,list,list,u,E,N,a,l,j,n,z1,z2, B,Nr,Rmax,m).Green_function;
  m.wave_function=non_local_wavefunction_matrix(original,list,list,u,E,N,a,l,j,n,z1,z2, B,Nr,Rmax,m).wave_function;
 for( int j=0; j<m.wave_function.size(); j++) {
            std::cout << "WAVE"<<" "<<m.wave_function[j].radius<<" "<<m.wave_function[j].real<<" "<<m.wave_function[j].imag << "  "<<endl;
        }
     


  for (int i=0; i<m.wave_function.size();i++)
 {
 
 outdata<< m.wave_function[i].radius<<"  "<<m.wave_function[i].real<<"  "<< m.wave_function[i].imag<<endl;
 }
 outdata.close();
  return m;
  
}




void printArray(double arr[], int size) {
    for ( int i = 0; i < size; i++ ) {
        cout << arr[i] << ' ';
    }
    cout << endl;
}


void print_matrix( double x[3][3] ) {

    for(unsigned int i=0; i< 3; i++) {
        for(unsigned int j=0; j<3; j++) {
            std::cout << x[i][j] << "  ";
        }
        std::cout << std::endl;
    }
}
int main()
{ 
 std::vector <double> u;
   int s;
   std::complex < double >** matrixx;
   //double**matrix_real;
   //double**matrix_imag;
   std::complex < double >** potential_matrix;
   std::complex < double >** potential_matrix_1;
 
 int start_s=clock(); 
 Final_return b;
double mu;
mu=40.0/41.0;
 non_local_wavefunction_matrixx(mu,5,60,40,1,1.5,20,20,0,0,0.1,30,b);
 
 

 int stop_s=clock();
 cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)<< "s"<<endl;
 
 //gsl_sf_hyperg_U (0.1, 0.1, 0.1);
 return 0;
 // matrixx=k_matrix(3,3);
 
 // printMatrix(matrixx,3);


}

