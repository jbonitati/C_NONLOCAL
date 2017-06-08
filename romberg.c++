 
// Integration with the romberg method :
// ------------------------------------- 

template <typename TYPE>   
TYPE romberg (const int N,const class array<TYPE> &F,const double x_debut,const double x_end)
{
  const int n = static_cast<int> (rint (log (N)/log (2)));
  const double h = x_end - x_debut;
  class array<TYPE> T(n+1,n+1);

  // T : table with the more and more precise value of the integral.
  // h : length of the interval of integration.

  int two_to_the_power_i = 1;
  T(0,0) = (F(0) + F(N))*h/2.0;
  for (int i = 1 ; i <= n ; i++)
  {
    two_to_the_power_i *= 2; 
    T(i,0) = 0.0;

    for (int j = 1 ; j < two_to_the_power_i ; j += 2)
    {
      int k = j*N/two_to_the_power_i;
      T(i,0) += F(k);
    }    
    T(i,0) = T(i-1,0)/2.0 + h*T(i,0)/two_to_the_power_i;
  }

  int four_to_the_power_i = 1;  
  for (int i = 1 ; i <= n ; i++)
  {
    four_to_the_power_i *= 4;
    int four_to_the_power_i_minus_one = four_to_the_power_i - 1;

    for (int j = 0 ; j <= n - i ; j++)
      T(j,i) = (four_to_the_power_i*T(j+1,i-1) - T(j,i-1))/four_to_the_power_i_minus_one;
  }

  return T(0,n);
}
