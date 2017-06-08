/**************************************
 * Legendre Polynomial
 * 
 * Kai Zhang, Duke University, 2011
 *************************************/
#include <iostream>
#include <fenv.h>   // enable floating point trap
#include <fstream>

using namespace std;
static void __attribute__((constructor)) trapfpe () {
 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}  // program will stop if an invalid number is generated

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

int main()
{
 int i,N;
 double x;
 double dx; // x increment

 
 cout << "Input order (integer n >= 0) of Legendre Polynomial:" << endl;
 cin >> N;
 cout << "Input x (double [-1, 1]):" << endl;
 cin >> x;
 cout << "Evaluate Legendre polynomial Pn(x) of order " << N << " at " << x << endl;

 cout << "P" << N << "(x=" << x << ")" << " = " << Legendre(N,x) << endl;

 ofstream os;
 os.open("px.txt");
 x=-1.0;
 dx = 0.01;
 while(x>=-1.0 && x<=1.0000000000001)
 {
  os << "x = " << x << "\t";
  for(i=0;i<=N;i++)
   os << "P" << i << "(x) = " << Legendre(i,x) << "\t";

  os << endl;
 
  x += dx;
 }
 os.close();

 return 0;
}
