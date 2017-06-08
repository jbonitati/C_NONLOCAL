//  21/09/2005 il y a inverse_iteration 


namespace Householder_QL
{
  template <typename TYPE> 
  void all_eigenpairs (class matrix<TYPE> &A,
		       class array<TYPE> &eigenvalues,
		       class vector_class<TYPE> **&eigenvectors);


  
  template <typename TYPE> 
  void all_eigenpairs_inverse_iteration (class matrix<TYPE> &A,
					 class array<TYPE> &eigenvalues,class vector_class<TYPE> **&eigenvectors,
					 const double prec);

  
  template <typename TYPE> 
  void one_eigenpairs_inverse_iteration (class matrix<TYPE> &A,
					 TYPE &eigenvalue_guess ,class vector_class<TYPE> &eigenvect_guess,double precision_eigenvect
					 ,double precision_eigenvalue);



 template <typename TYPE> 
  void all_eigenpairs2 (class matrix<TYPE> &A,
		       class array<TYPE> &eigenvalues,
		       class vector_class<TYPE> **&eigenvectors);

 template <typename TYPE> 
  void all_eigenpairs3 (class matrix<TYPE> &A,
		       class array<TYPE> &eigenvalues,
		       class vector_class<TYPE> **&eigenvectors);


 template <typename TYPE> 
 void all_eigenpairs3_check (class matrix<TYPE> &A,
			     class array<TYPE> &eigenvalues,
			     class vector_class<TYPE> **&eigenvectors,class matrix<TYPE> &B);



template <typename TYPE>
  void all_eigenpairs4 (class matrix<TYPE> &A,
                       class array<TYPE> &eigenvalues,
                       class vector_class<TYPE> **&eigenvectors);
                                                                                                                                             
                                                                                                                                             



  template <typename TYPE> 
  void all_eigenpairs (const int N,
		       const class array<TYPE> &diag,
		       const class array<TYPE> &off_diag,
		       class array<TYPE> &eigenvalues,
		       class vector_class<TYPE> **&eigenvectors);
  


  template <typename TYPE> 
  void all_eigenvalues (class matrix<TYPE> &A,
			class array<TYPE> &eigenvalues);

  template <typename TYPE> 
  void all_eigenvalues (const int N,
			const class array<TYPE> &diag,
			const class array<TYPE> &off_diag,
			class array<TYPE> &eigenvalues);

  template <typename TYPE> 
  void eigen_sort (const int low,const int high,const int N,const bool with_vector,
		   class array<TYPE> &eigenvalues,
		   class vector_class<TYPE> * eigenvectors[]);


 template <typename TYPE> 
  void eigen_sort2 (const int low,const int high,const int N,const bool with_vector,
		   class array<TYPE> &eigenvalues,
		   class vector_class<TYPE> * eigenvectors[]);

template <typename TYPE> 
  void eigen_sort3 (const int low,const int high,const int N,const bool with_vector,
		   class array<TYPE> &eigenvalues,
		   class vector_class<TYPE> * eigenvectors[]);



template <typename TYPE>
  void eigen_sort4 (const int low,const int high,const int N,const bool with_vector,
                   class array<TYPE> &eigenvalues,
                   class vector_class<TYPE> * eigenvectors[]);
                                                                                                                                             



  template <typename TYPE> 
  void tridiag (class matrix<TYPE> &A,
		class array<TYPE> &diag,
		class array<TYPE> &off_diag);

  template <typename TYPE> 
  void QL_diag (const bool with_vector,const int dimension,
		class matrix<TYPE> *const A_ptr,
		class array<TYPE> &diag,
		class array<TYPE> &off_diag);

 
  template <typename TYPE> 
  void QL_diag_check (const bool with_vector,const int dimension,
		      class matrix<TYPE> *const A_ptr,
		      class array<TYPE> &diag,
		      class array<TYPE> &off_diag,class matrix<TYPE> & B);

}




// Initial values :
// low = 0.
// high = N - 1.

template <typename TYPE> 
void Householder_QL::eigen_sort (const int low,const int high,const int N,const bool with_vector,
				 class array<TYPE> &eigenvalues,
				 class vector_class<TYPE> * eigenvectors[])
{
  int i_sort = low;
  int j_sort = high;
  const double pivot_eigenvalue = real (eigenvalues((i_sort + j_sort)/2));

  do
  {
    while ((real (eigenvalues(i_sort)) < pivot_eigenvalue) )
      i_sort++;

    while ((real (eigenvalues(j_sort)) > pivot_eigenvalue) )
      j_sort--;
     
    if (i_sort <= j_sort)
    {
      const TYPE temp_eigenvalue = eigenvalues(i_sort);
      eigenvalues(i_sort) = eigenvalues(j_sort);
      eigenvalues(j_sort) = temp_eigenvalue;

      if (with_vector == true)
      {
	class vector_class<TYPE> *const temp_eig_vec = eigenvectors[i_sort];
	eigenvectors[i_sort] = eigenvectors[j_sort];
	eigenvectors[j_sort] = temp_eig_vec;
      }

      i_sort++;
      j_sort--;

      
    }
  }
  while (i_sort <= j_sort);
  
  if (low < j_sort) eigen_sort (low,j_sort,N,with_vector,eigenvalues,eigenvectors);

  if (i_sort < high) eigen_sort (i_sort,high,N,with_vector,eigenvalues,eigenvectors); 
}


template <typename TYPE> 
void Householder_QL::eigen_sort2 (const int low,const int high,const int N
				  ,const bool with_vector,class array<TYPE> &eigenvalues
				  ,class vector_class<TYPE> * eigenvectors[])
{
  int i_sort = low;
  int j_sort = high;
  const double pivot_eigenvalue = abs (real (eigenvalues((i_sort + j_sort)/2)));

  do
  {
    while ( ( abs( real (eigenvalues(i_sort))  ) < pivot_eigenvalue) )
      i_sort++;

    while ( ( abs(real (eigenvalues(j_sort))  ) > pivot_eigenvalue) )
      j_sort--;
     
    if (i_sort <= j_sort)
    {
      const TYPE temp_eigenvalue = eigenvalues(i_sort);
      eigenvalues(i_sort) = eigenvalues(j_sort);
      eigenvalues(j_sort) = temp_eigenvalue;

      if (with_vector == true)
      {
	class vector_class<TYPE> *const temp_eig_vec = eigenvectors[i_sort];
	eigenvectors[i_sort] = eigenvectors[j_sort];
	eigenvectors[j_sort] = temp_eig_vec;
      }

      i_sort++;
      j_sort--;

    }
  }
  while (i_sort <= j_sort);
  
  if (low < j_sort) eigen_sort2 (low,j_sort,N,with_vector,eigenvalues,eigenvectors);

  if (i_sort < high) eigen_sort2 (i_sort,high,N,with_vector,eigenvalues,eigenvectors); 
}


template <typename TYPE> 
void Householder_QL::eigen_sort3 (const int low,const int high,const int N                         
				  ,const bool with_vector,class array<TYPE> &eigenvalues
				  ,class vector_class<TYPE> * eigenvectors[])
{
  int i_sort = low;
  int j_sort = high;
  const double pivot_eigenvalue = abs ( (eigenvalues((i_sort + j_sort)/2)));

  do
  {
    while (  abs((eigenvalues(i_sort))) < pivot_eigenvalue )
      {
	i_sort++;
      }


    while (  abs((eigenvalues(j_sort))) > pivot_eigenvalue )
      {
	j_sort--;
      }
      

    if (i_sort <= j_sort)
    {
      const TYPE temp_eigenvalue = eigenvalues(i_sort);
      eigenvalues(i_sort) = eigenvalues(j_sort);
      eigenvalues(j_sort) = temp_eigenvalue;

      if (with_vector == true)
      {
	class vector_class<TYPE> *const temp_eig_vec = eigenvectors[i_sort];
	eigenvectors[i_sort] = eigenvectors[j_sort];
	eigenvectors[j_sort] = temp_eig_vec;
      }

      i_sort++;
      j_sort--;


    }
  }
  while (i_sort <= j_sort);
  
  if (low < j_sort) eigen_sort3 (low,j_sort,N,with_vector,eigenvalues,eigenvectors);

  if (i_sort < high) eigen_sort3 (i_sort,high,N,with_vector,eigenvalues,eigenvectors); 
}


template <typename TYPE>
void Householder_QL::eigen_sort4 (const int low,const int high,const int N
                                  ,const bool with_vector,class array<TYPE> &eigenvalues
                                  ,class vector_class<TYPE> * eigenvectors[])
{
  int i_sort = low;
  int j_sort = high;


  const double pivot_eigenvalue = -abs ( (eigenvalues((i_sort + j_sort)/2)));

  do
  {

    while ( ( -abs((eigenvalues(i_sort))) < pivot_eigenvalue ))
	    
      {
	i_sort++;
      }

                                                                                                                                         
    while ( ( -abs( (eigenvalues(j_sort))) > pivot_eigenvalue  ))
	   
      {j_sort--;}

  

  if (i_sort <= j_sort)
    {
      const TYPE temp_eigenvalue = eigenvalues(i_sort);
      eigenvalues(i_sort) = eigenvalues(j_sort);
      eigenvalues(j_sort) = temp_eigenvalue;
                                                                                                                                             
      if (with_vector == true)
      {
        class vector_class<TYPE> *const temp_eig_vec = eigenvectors[i_sort];
        eigenvectors[i_sort] = eigenvectors[j_sort];
        eigenvectors[j_sort] = temp_eig_vec;
      }
                                                                                                                                             
      i_sort++;
      j_sort--;
                                                                                                                                             
    }
  }

  while (i_sort <= j_sort);
                                                                                                                                             
  if (low < j_sort) eigen_sort4 (low,j_sort,N,with_vector,eigenvalues,eigenvectors);
                                                                                                                                             
  if (i_sort < high) eigen_sort4 (i_sort,high,N,with_vector,eigenvalues,eigenvectors);
  }
                                                                                                                                             



template <typename TYPE> 
void Householder_QL::tridiag (class matrix<TYPE> &A,
			      class array<TYPE> &diag,
			      class array<TYPE> &off_diag)
{
  for (int i = A.dim-1 ; i >= 1 ; i--) 
  {
    int l = i - 1;
    double scale = 0.0;
    TYPE h = 0.0;
    
    if (l > 0) 
    {
      for (int k = 0 ; k <= l ; k++)
	scale += abs (A(i,k));
      
      if (scale == 0.0)
	off_diag(i) = A(i,l);
      else
      {
	for (int k = 0 ; k <= l ; k++) 
	{
	  A(i,k) /= scale;
	  h += A(i,k)*A(i,k);
	}
	
	TYPE f = A(i,l);
	TYPE g = (real (f) > 0) ? (-sqrt (h)) : (sqrt (h));
	off_diag(i) = scale*g;
	h -= f*g;
	A(i,l) = f-g;
	f = 0.0;
	for (int j = 0 ; j <= l ; j++)
	{                       
	  A(j,i) = A(i,j)/h;
	  g = 0.0;
	  
	  for (int k = 0 ; k <= j ; k++)
	    g += A(j,k)*A(i,k);
	  for (int k = j+1 ; k <= l ; k++)
	    g += A(k,j)*A(i,k);
	  
	  off_diag(j) = g/h;
	  f += off_diag(j)*A(i,j);
	}
	
	TYPE hh = f/(h+h);
	for (int j = 0 ; j <= l ; j++)
	{
	  f = A(i,j);
	  off_diag(j) = g = off_diag(j) - hh*f;
	  for (int k = 0 ; k <= j ; k++)
	    A(j,k) -= f*off_diag(k) + g*A(i,k);
	}
      }
    }
    else
      off_diag(i) = A(i,l);

    diag(i) = h;
  }

  diag(0) = off_diag(0) = 0.0;

  for (int i = 0 ; i < A.dim ; i++)
  {
    int l = i-1;
    if (diag(i) != 0.0) 
      for (int j = 0 ; j <= l ; j++)
      {
	TYPE g = 0.0;
	for (int k = 0 ; k <= l ; k++)
	  g += A(i,k)*A(k,j);
	for (int k = 0 ; k <= l ; k++)
	  A(k,j) -= g*A(k,i);
      }
    
    diag(i) = A(i,i);
    A(i,i) = 1.0;

    for (int j = 0 ; j <= l ; j++)
      A(j,i) = A(i,j) = 0.0;
  }
}






template <typename TYPE> 
void Householder_QL::QL_diag (const bool with_vector,const int dimension,class matrix<TYPE> *const A_ptr,
			      class array<TYPE> &diag,class array<TYPE> &off_diag)
{   
  for (int l = 0 ; l < dimension ; l++) 
  {
    int iter = 0,m; 
    do 
    {
      for (m = l ; m < dimension - 1 ; m++) 
      {
	double dd = abs (diag(m)) + abs (diag(m+1)); 
       	if ( abs (off_diag(m))<1e-10) 
	  //  if ((abs (off_diag(m)) + dd) == dd)  

	  break; 
      }
      
      if (m != l) 
      {
	if (iter++ == 1000) 
	{
	  cout<<"Too many iterations."<<endl;
          cout<<"|diag(m)|+|diag(m+1)|:"<<abs (diag(m)) + abs (diag(m+1))<<endl;
          cout<<"|off_diag(m)|:"<<abs (off_diag(m))<<endl;


	  // for (int xc=0;xc<=dimension-1;xc++)
	  //   for (int xc2=xc;xc2<=dimension-1;xc2++)
	  //     {
	  // 	cout<<xc<<" "<<xc2<<" "<<(*A_ptr)(xc,xc2)<<endl;
	  //     }

	  exit (1);
	} 
	
	TYPE g = (diag(l+1) - diag(l))/(2.0*off_diag(l));
	TYPE r = sqrt (1.0 + g*g);
	
	(abs (g + r) > abs (g - r))
	  ? (g = diag(m) - diag(l) + off_diag(l)/(g + r)) 
	  : (g = diag(m) - diag(l) + off_diag(l)/(g - r)); 
	
	TYPE s = 1.0, c = 1.0, p = 0.0;
	
	int i;
	for (i = m-1 ; i >= l ; i--) 
	{
	  TYPE f = s*off_diag(i),b = c*off_diag(i);
	  off_diag(i+1) = (r = sqrt (f*f + g*g)); 
	  
	  if (r == 0.0) 
	  {
	    diag(i+1) -= p; 
	    off_diag(m) = 0.0; 
	    break; 
	  }
	  
	  s = f/r; 
	  c = g/r; 
	  g = diag(i+1) - p; 
	  r = (diag(i) - g)*s + 2.0*c*b; 
	  diag(i+1) = g + (p = s*r); 
	  g = c*r - b;

	  if (with_vector == true)
	  {
	    class matrix<TYPE> &A = *A_ptr;
	    for (int k = 0 ; k < dimension ; k++)
	    {
	      f = A(k,i+1);
	      A(k,i+1) = s*A(k,i) + c*f;
	      A(k,i) = c*A(k,i) - s*f;
	    }
	  } 
	}

	
	if ((r == 0.0) && (i >= l)) 
	  continue;
 
	diag(l) -=  p ; 
	off_diag(l) = g ; 
	off_diag(m) = 0.0 ; 	
      } 
    }
    while (m != l); 
  }
}  











template <typename TYPE> 
void Householder_QL::QL_diag_check (const bool with_vector,const int dimension,class matrix<TYPE> *const A_ptr,
				    class array<TYPE> &diag,class array<TYPE> &off_diag,class matrix<TYPE> & B)
{   
  for (int l = 0 ; l < dimension ; l++) 
  {
    int iter = 0,m; 
    do 
    {
      for (m = l ; m < dimension - 1 ; m++) 
      {
	double dd = abs (diag(m)) + abs (diag(m+1)); 
       	if ( abs (off_diag(m))<1e-10) 
	  //  if ((abs (off_diag(m)) + dd) == dd)  

	  break; 
      }
      
      if (m != l) 
      {
	if (iter++ == 1300) 
	{
	  cout<<"HERE Too many iterations."<<endl;
         

	  ofstream out_file ("matrix_fail");
          cout<<"matrix fail "<<endl;
	  out_file.precision(15);

	  for (int zw = 0 ; zw < dimension ; zw++) 
	    for (int as = 0 ; as < dimension ; as++) 
	      {
		out_file<<B(zw,as)<<endl;	      
	      }

	  out_file.close();

	  for (int zw = 0 ; zw < dimension ; zw++)
            for (int as = 0 ; as < dimension ; as++)
              {
		cout<<zw<<" "<<as<<" "<<B(zw,as)<<endl;
	      }


	  exit (1);
	} 
	
	TYPE g = (diag(l+1) - diag(l))/(2.0*off_diag(l));
	TYPE r = sqrt (1.0 + g*g);
	
	(abs (g + r) > abs (g - r))
	  ? (g = diag(m) - diag(l) + off_diag(l)/(g + r)) 
	  : (g = diag(m) - diag(l) + off_diag(l)/(g - r)); 
	
	TYPE s = 1.0, c = 1.0, p = 0.0;
	
	int i;
	for (i = m-1 ; i >= l ; i--) 
	{
	  TYPE f = s*off_diag(i),b = c*off_diag(i);
	  off_diag(i+1) = (r = sqrt (f*f + g*g)); 
	  
	  if (r == 0.0) 
	  {
	    diag(i+1) -= p; 
	    off_diag(m) = 0.0; 
	    break; 
	  }
	  
	  s = f/r; 
	  c = g/r; 
	  g = diag(i+1) - p; 
	  r = (diag(i) - g)*s + 2.0*c*b; 
	  diag(i+1) = g + (p = s*r); 
	  g = c*r - b;

	  if (with_vector == true)
	  {
	    class matrix<TYPE> &A = *A_ptr;
	    for (int k = 0 ; k < dimension ; k++)
	    {
	      f = A(k,i+1);
	      A(k,i+1) = s*A(k,i) + c*f;
	      A(k,i) = c*A(k,i) - s*f;
	    }
	  } 
	}

	
	if ((r == 0.0) && (i >= l)) 
	  continue;
 
	diag(l) -=  p ; 
	off_diag(l) = g ; 
	off_diag(m) = 0.0 ; 	
      } 
    }
    while (m != l); 
  }
}  





// A is changed to N vectors then destroyed. eigenvectors must not be allocated.
// eigenpairs ordered.
template <typename TYPE> 
void Householder_QL::all_eigenpairs (class matrix<TYPE> &A,
				     class array<TYPE> &eigenvalues,
				     class vector_class<TYPE> **&eigenvectors)
{ 
  class array<TYPE> off_diag(A.dim);

  tridiag (A,eigenvalues,off_diag);

  for (unsigned int i = 1 ; i < A.dim ; i++) off_diag(i-1) = off_diag(i);
  off_diag(A.dim-1) = 0.0;

  QL_diag (true,A.dim,&A,eigenvalues,off_diag);
  A.matrix_to_vectors (eigenvectors);

  eigen_sort (0,A.dim-1,A.dim,true,eigenvalues,eigenvectors);
}




template <typename TYPE> 
void Householder_QL::all_eigenpairs_inverse_iteration (class matrix<TYPE> &A,
						       class array<TYPE> &eigenvalues,class vector_class<TYPE> **&eigenvectors,
						       const double prec)
{ 
  const double A_inf = A.infinite_norm ();
  A /= A_inf;

  const class matrix<TYPE> AA = A;

  const unsigned int N = A.dim;

  for (unsigned int i = 0 ; i < N ; i++) 
    for (unsigned int j = 0 ; j < i  ; j++) 
      if (inf_norm (A(i,j)) < prec) A(i,j) = A(j,i) = 0.0;

  all_eigenpairs (A,eigenvalues,eigenvectors);

  for (unsigned int i = 0 ; i < N ; i++)
    for (unsigned int loop = 0 ; loop < 2 ; loop++)
    {
      class matrix<TYPE> A_eig_i = AA;
      for (unsigned int ii = 0 ; ii < N ; ii++) { A_eig_i(ii,ii) -= eigenvalues(i);}

      (*(eigenvectors[i]))=linear_system_solution (A_eig_i,(*(eigenvectors[i])));
      (*(eigenvectors[i])).normalization();

      eigenvalues(i) = (*(eigenvectors[i])) *(AA*(*(eigenvectors[i])));


     //  A.eigenvector(i) = linear_system_solution (A_eig_i,A.eigenvector(i));
//       A.eigenvector(i).normalization ();
//       eigenvalues(i) = A.eigenvector(i)*(AA*A.eigenvector(i));

    }

  for (unsigned int i = 0 ; i < N ; i++) eigenvalues(i) *= A_inf;

 //  for (unsigned int i = 0 ; i < N ; i++)
//     {
//       cout<<"val pro "<<eigenvalues(i)<<endl;
      
//     }


}


template <typename TYPE> 
void Householder_QL::one_eigenpairs_inverse_iteration (class matrix<TYPE> &A,
						       TYPE &eigenvalue_guess ,class vector_class<TYPE> &eigenvect_guess,double precision_eigenvect
						       ,double precision_eigenvalue)
{

  const unsigned int N = A.dim;

  class vector_class<TYPE> previous_vect(N);
  TYPE previous_eigenvalue;
  
  class vector_class<TYPE> residu_vect(N);
  TYPE residu_eigenvalue;

  do
    {
      previous_vect=eigenvect_guess;
      previous_eigenvalue=eigenvalue_guess;
     
      
      class matrix<TYPE> A_eig_i = A;
      
      for (unsigned int ii = 0 ; ii < N ; ii++)
	{
	  A_eig_i(ii,ii) -= eigenvalue_guess;
	}
      
      class vector_class<TYPE> temp_vect_guess=linear_system_solution (A_eig_i,eigenvect_guess);
	  
      temp_vect_guess.normalization();
	  
      TYPE temp_eigenvalue=eigenvect_guess*(A*eigenvect_guess);
      
      eigenvect_guess=temp_vect_guess;
      eigenvalue_guess=temp_eigenvalue;
      
      if ( real(eigenvect_guess.coef_max())*real(previous_vect.coef_max())<0.0) {eigenvect_guess*=-1.0;}
      
      residu_vect= previous_vect-eigenvect_guess;
      residu_eigenvalue= previous_eigenvalue-eigenvalue_guess;

      cout<<" residu_eigenvect et residu_eigenvalue "<<residu_vect.dot_product_norm()<<" "<<residu_eigenvalue<<endl;
      
    }

  while (residu_vect.infinite_norm()>precision_eigenvect  || abs (residu_eigenvalue)>precision_eigenvalue );

  cout<<"final value "<<eigenvalue_guess<<endl;
  cout<<" "<<endl;

}






template <typename TYPE> 
void Householder_QL::all_eigenpairs2 (class matrix<TYPE> &A,
				     class array<TYPE> &eigenvalues,
				     class vector_class<TYPE> **&eigenvectors)
{ 
  class array<TYPE> off_diag(A.dim);

  tridiag (A,eigenvalues,off_diag);

  for (unsigned int i = 1 ; i < A.dim ; i++) off_diag(i-1) = off_diag(i);
  off_diag(A.dim-1) = 0.0;

  QL_diag (true,A.dim,&A,eigenvalues,off_diag);
  A.matrix_to_vectors (eigenvectors);

  eigen_sort2 (0,A.dim-1,A.dim,true,eigenvalues,eigenvectors);
}



template <typename TYPE> 
void Householder_QL::all_eigenpairs3 (class matrix<TYPE> &A,
				     class array<TYPE> &eigenvalues,
				     class vector_class<TYPE> **&eigenvectors)
{ 
  class array<TYPE> off_diag(A.dim);

  tridiag (A,eigenvalues,off_diag);

  for (unsigned int i = 1 ; i < A.dim ; i++) off_diag(i-1) = off_diag(i);
  off_diag(A.dim-1) = 0.0;

  QL_diag (true,A.dim,&A,eigenvalues,off_diag);
  A.matrix_to_vectors (eigenvectors);

  eigen_sort3 (0,A.dim-1,A.dim,true,eigenvalues,eigenvectors);
}



template <typename TYPE> 
void Householder_QL::all_eigenpairs3_check (class matrix<TYPE> &A,
					    class array<TYPE> &eigenvalues,
					    class vector_class<TYPE> **&eigenvectors,class matrix<TYPE> &B)
{ 
  class array<TYPE> off_diag(A.dim);

  tridiag (A,eigenvalues,off_diag);

  for (unsigned int i = 1 ; i < A.dim ; i++) off_diag(i-1) = off_diag(i);
  off_diag(A.dim-1) = 0.0;

  QL_diag_check (true,A.dim,&A,eigenvalues,off_diag,B);
  A.matrix_to_vectors (eigenvectors);

  eigen_sort3 (0,A.dim-1,A.dim,true,eigenvalues,eigenvectors);
}


                                                                                                                                             
template <typename TYPE>
void Householder_QL::all_eigenpairs4 (class matrix<TYPE> &A,
                                     class array<TYPE> &eigenvalues,
                                     class vector_class<TYPE> **&eigenvectors)
{
  class array<TYPE> off_diag(A.dim);
                                                                                                                                             
  tridiag (A,eigenvalues,off_diag);
                                                                                                                                             
  for (unsigned int i = 1 ; i < A.dim ; i++) off_diag(i-1) = off_diag(i);
  off_diag(A.dim-1) = 0.0;
                                                                                                                                             
  QL_diag (true,A.dim,&A,eigenvalues,off_diag);
  A.matrix_to_vectors (eigenvectors);
                                                                                                                                             
  eigen_sort4 (0,A.dim-1,A.dim,true,eigenvalues,eigenvectors);
}



// eigenvectors must not be allocated. Tridiagonal matrices only.
// eigenpairs ordered.
template <typename TYPE> 
void Householder_QL::all_eigenpairs (const int N,
				     const class array<TYPE> &diag,
				     const class array<TYPE> &off_diag,
				     class array<TYPE> &eigenvalues,
				     class vector_class<TYPE> **&eigenvectors)
{ 
  class array<TYPE> off_diag_temp (N);

  for (unsigned int i = 0 ; i < N ; i++) off_diag_temp(i) = off_diag(i);
  for (unsigned int i = 0 ; i < N ; i++) eigenvalues(i) = diag(i);
 
  class matrix<TYPE> Q(N);
  Q.identity ();
  QL_diag (true,N,&Q,eigenvalues,off_diag_temp);
  Q.matrix_to_vectors (eigenvectors);

  eigen_sort (0,N-1,N,true,eigenvalues,eigenvectors);
}




// No eigenvectors calculated. Only eigenvalues.
// eigenvalues ordered.
template <typename TYPE> 
void Householder_QL::all_eigenvalues (class matrix<TYPE> &A,
				      class array<TYPE> &eigenvalues)
{
  class array<TYPE> off_diag(A.dim);

  tridiag (A,eigenvalues,off_diag);

  for (unsigned int i = 1 ; i < A.dim ; i++) off_diag(i-1) = off_diag(i);
  off_diag(A.dim-1) = 0.0;

  QL_diag (false,A.dim,&A,eigenvalues,off_diag);

  class vector_class<TYPE> * *const NOTHING0 = 0;
  eigen_sort (0,A.dim - 1,A.dim,false,eigenvalues,NOTHING0);
}



// No eigenvectors calculated. Only eigenvalues. Tridiagonal matrices only.
// eigenvalues ordered.
template <typename TYPE> 
void Householder_QL::all_eigenvalues (const int N,
				      const class array<TYPE> &diag,
				      const class array<TYPE> &off_diag,
				      class array<TYPE> &eigenvalues)
{
  class array<TYPE> off_diag_temp (N);

  for (unsigned int i = 0 ; i < N ; i++) off_diag_temp(i) = off_diag(i);
  
  class matrix<TYPE> *const NOTHING0 = 0;					    

  for (unsigned int i = 0 ; i < N ; i++) eigenvalues(i) = diag(i);

  QL_diag (false,N,NOTHING0,eigenvalues,off_diag_temp);

  class vector_class<TYPE> * *const NOTHING1 = 0;
  eigen_sort (0,N-1,N,false,eigenvalues,NOTHING1);
}
