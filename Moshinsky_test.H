bool is_9j_trivial_zero (const double  j1,const double  j2,const double J12,
                                           const double  j3,const double  j4,const double J34,
                                           const double J13,const double J24,const double J) 
{
  if (is_it_triangle (j1,j2,J12) == false) return true;
  if (is_it_triangle (j3,j4,J34) == false) return true;
  if (is_it_triangle (J13,J24,J) == false) return true; 

  if (is_it_triangle (j1,j3,J13) == false) return true;
  if (is_it_triangle (j2,j4,J24) == false) return true;
  if (is_it_triangle (J12,J34,J) == false) return true; 

  if (is_it_integer (j1 + j2 + J12) == false) return true;
  if (is_it_integer (j3 + j4 + J34) == false) return true;
  if (is_it_integer (J13 + J24 + J) == false) return true;
  
  if (is_it_integer (j1 + j3 + J13) == false) return true;
  if (is_it_integer (j2 + j4 + J24) == false) return true;
  if (is_it_integer (J12 + J34 + J) == false) return true;

  return false;
}


double ln_n_fact(int n) // returns ln (n!)
{
  double res=0.0;

  for (int i=2;i<=n;i++)
    {
      res+=log(i*1.0);
    }

  return res;

}




double func_mosh_term1(int n1,int n2,int k,int s1,int s2,int s3,int q1,int q2,int q3)
{
  double ln_res=0;

  ln_res=//ln_n_fact(n1)+ln_n_fact(n2)+
    k*1.0*log(2.0)
    -ln_n_fact(s1)-ln_n_fact(s2)-ln_n_fact(s3)-ln_n_fact(q1)-ln_n_fact(q2)-ln_n_fact(q3);

  double res=exp(ln_res);

  return res;
}

double func_G(int J,int k)
{
  // double G = (2*J+1)*pow(2.0,J)*fact_k*factorial (0.5*(k + J))/factorial (k + J + 1)/factorial (0.5*(k - J));

  int half1=0.5*(k + J);
  int half2=0.5*(k - J);
  

  double ln_G=log(2*J+1)+J*log(2.0)+ln_n_fact(k)+ln_n_fact(half1)-ln_n_fact(k+J+1)-ln_n_fact(half2);

  double G=exp(ln_G);

  return G;

}


double Moshinsky (const int n1,const int l1,const int n2,const int l2,const int N,const int L,const int n,const int l,const int lambda)
{
  const double fact_n1 = factorial(n1), fact_n2 = factorial(n2), fact_n1_x_fact_n2 = fact_n1*fact_n2;
 
  double mosh = 0.0;  
  for (int F1 = 0 ; F1 <= l1 ; F1++)
    {
      const int f1 = l1 - F1;
      const double coef_bin_2l1_2F1 = factorial(2*l1)/factorial(2*F1)/factorial(2*(l1 - F1));

      for (int F2 = 0 ; F2 <= l2 ; F2++)
        {
          const int f2 = l2 - F2;
          const double coef_bin_2l2_2F2 = factorial(2*l2)/factorial(2*F2)/factorial(2*(l2 - F2)),coef_bin_prod = coef_bin_2l1_2F1*coef_bin_2l2_2F2;

          for (int s1 = 0 ; s1 <= n1 ; s1++)
            for (int s2 = 0 ; s2 <= n1 ; s2++)
              {
                const int s3 = n1 - s1 - s2;
                if (s3 >= 0)
                  {
                    const double fact_s1=factorial(s1),fact_s2=factorial(s2),fact_s3=factorial(s3),fact_s1_x_fact_s2_x_fact_s3=fact_s1*fact_s2*fact_s3;

                    for (int q1 = 0 ; q1 <= n2 ; q1++)
                      for (int q2 = 0 ; q2 <= n2 ; q2++)
                        {
                          const int q3 = n2 - q1 - q2;
                          if (q3 >= 0)
                            {
                              const int k = s3 + q3;
                              const double fact_k = factorial (k),fact_n1_x_fact_n2_x_two_pow_k = fact_n1_x_fact_n2*pow(2.0,k),fact_q1 = factorial(q1);
                              const double fact_q2=factorial(q2),fact_q3=factorial(q3),fact_q1_x_fact_q2_x_fact_q3=fact_q1*fact_q2*fact_q3;
                              
                              if ((F1 + F2 + 2*s1 + 2*q1 + s3 + q3 == 2*N + L) && (f1 + f2 + 2*s2 + 2*q2 + s3 + q3 == 2*n + l))
                                for (int L1 = abs (F1 - F2) ; L1 <= F1 + F2 ; L1++)
                                  {
                                    const double CG1 = Clebsch_Gordan (F1,0,F2,0,L1,0);
                                    
                                    for (int L2 = abs (f1 - f2) ; L2 <= f1 + f2 ; L2++)
                                      {
                                        const double CG2 = Clebsch_Gordan (f1,0,f2,0,L2,0),Wig1_9j = Wigner_9j(F1,f1,l1,F2,f2,l2,L1,L2,lambda);
                                        const double CG1_x_CG2_x_Wig1_9j = CG1*CG2*Wig1_9j;

                                        for (int J = k ; J >= 0; J -= 2)
                                          {
                                            const int mosh_phase = ((s3 + f2) % 2 == 0) ? (1) : (-1);
					    //    const double mosh_term1 = fact_n1_x_fact_n2_x_two_pow_k/fact_s1_x_fact_s2_x_fact_s3/fact_q1_x_fact_q2_x_fact_q3;
					     const double mosh_term1 = pow(2.0,k)/fact_s1_x_fact_s2_x_fact_s3/fact_q1_x_fact_q2_x_fact_q3;

					    double mosh_term1_bis=func_mosh_term1(n1,n2,k,s1,s2,s3,q1,q2,q3);

					
						

                                            const double mosh_term2 = sqrt(coef_bin_prod*(2*L1 + 1)*(2*L2 + 1)*(2*J + 1)*(2*lambda + 1));
                                            const double CG3 = Clebsch_Gordan (J,0,L1,0,L,0),CG4 = Clebsch_Gordan (J,0,L2,0,l,0);
                                            const double G = (2*J+1)*pow(2.0,J)*fact_k*factorial (0.5*(k + J))/factorial (k + J + 1)/factorial (0.5*(k - J));
                                            const double mosh_term3 = pow(0.5,f1 + 2*s2 + s3 + f2 + 2*q2 + q3),Wig2_9j = Wigner_9j(J,J,0,L1,L2,lambda,L,l,lambda); 
                                         


					    if ( triang(J,L2,l)==true   && triang(J,L1,L)==true &&  triang(F1,F2,L1)==true && triang(f1,f2,L2)
						 && is_9j_trivial_zero(J,J,0,L1,L2,lambda,L,l,lambda)==false && is_9j_trivial_zero (F1,f1,l1,F2,f2,l2,L1,L2,lambda)==false )
					      {
						mosh += mosh_phase*mosh_term1*mosh_term2*CG1_x_CG2_x_Wig1_9j*CG3*CG4*G*mosh_term3*Wig2_9j;
						
						if (abs(mosh_term1_bis-mosh_term1)>1e-5)
						  {
						    cout<<" ACCCCCCCCCCCCCCCCHTUNGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG  mosh_term1 "<<mosh_term1_bis<<" "<<mosh_term1<<endl;
						    cout<<"pb "<<n1<<" "<<n2<<" "<<N<<" "<<n<<endl;
						    cout<<"mosh tmp "<<mosh_phase*mosh_term1*mosh_term2*CG1_x_CG2_x_Wig1_9j*CG3*CG4*G*mosh_term3*Wig2_9j<<endl;
						    cout<<"G   "<<G<<" "<<func_G(J,k)<<endl;
						  }
					      }

                                          }}}}}}}}}
  
  const int phase = ((n1 + n2 + N + n) % 2 == 0) ? (1) : (-1);
  const double factor1 = pow(2.0,0.5*(n1 + n2 + l1 + l2 - N - n - L - l)),factor2 = pow (2.0,n - N + 0.5*(l - L));


  //  const double factor3 = sqrt(factorial(N)*factorial(n)*double_factorial(2*L + 2*N + 1)*double_factorial(2*l + 2*n + 1)
  //                          /fact_n1/fact_n2/double_factorial(2*l1 + 2*n1 + 1)/double_factorial(2*l2 + 2*n2 + 1));

  const double factor3 = sqrt(factorial(N)*factorial(n)*double_factorial(2*L + 2*N + 1)*double_factorial(2*l + 2*n + 1)
			      /double_factorial(2*l1 + 2*n1 + 1)/double_factorial(2*l2 + 2*n2 + 1));




  mosh *= phase*factor1*factor2*factor3*(2*l1 + 1)*(2*l2 + 1);

  cout<<"fact "<<phase*factor1*factor2*factor3*(2*l1 + 1)*(2*l2 + 1)<<endl;
  return mosh;
}
