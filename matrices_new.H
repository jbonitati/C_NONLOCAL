template <typename SCALAR_TYPE>
class matrix
{
public :
  matrix (const unsigned int dim1) : dim(dim1)
    {
      deleted = false;

      table = new SCALAR_TYPE * [dim];
      
      for (unsigned int i = 0 ; i < dim ; i++)
	{table[i] = new SCALAR_TYPE [dim];}

      row_permutation_private = new unsigned int [dim];

      is_it_LU_decomposed_private = false;


    }
  
  matrix (const class matrix &A) : dim (A.dim) ,is_it_LU_decomposed_private (A.is_it_LU_decomposed ())
    {
      
      deleted = false;

      row_permutation_private = new unsigned int [dim];

      for (unsigned int i = 0 ; i < dim ; i++) 
	row_permutation_private[i] = A.row_permutation(i);


      table = new SCALAR_TYPE * [dim];
      
      for (unsigned int i = 0 ; i < dim ; i++)
	table[i] = new SCALAR_TYPE [dim];

       for (unsigned int i = 0 ; i < dim ; i++) 
	 for (unsigned int j = 0 ; j < dim ; j++)
	   table[i][j] = A(i,j);
    }

  ~matrix (void)
    {
      if (deleted == false)
      {
	for (unsigned int i = 0 ; i < dim ; i++)
	  delete [] table[i];

	delete [] table;

	delete [] row_permutation_private;
      }
    }
  
  class matrix& operator += (const class matrix &);
  class matrix& operator -= (const class matrix &);
  class matrix& operator *= (const SCALAR_TYPE);
  class matrix& operator /= (const SCALAR_TYPE);

  const class matrix& operator = (const class matrix &);
  class matrix operator + (const class matrix &);
  class matrix operator - (const class matrix &);
  class matrix operator * (const class matrix &); 
  class matrix operator * (const SCALAR_TYPE) const; 
  class matrix operator / (const SCALAR_TYPE); 

  SCALAR_TYPE & operator() (const unsigned int,const unsigned int) const;
  void swap_vectors (const int unsigned i,const unsigned int j);
  

  void zero (void);
  void identity (void);
  void transpose (void);
  double infinite_norm (void) const;	  
  complex<double> trace (void) const;
  void matrix_to_vectors (class vector_class<SCALAR_TYPE> **&eigenvectors);

  void LU_decompose (void);
  
  const unsigned int dim;

  bool is_it_LU_decomposed (void) const
  {
    return is_it_LU_decomposed_private;
  }
  
  unsigned int & row_permutation (const unsigned int i) const
  {
    return row_permutation_private[i];
  }




private:
  bool deleted;
  SCALAR_TYPE **table;
  bool is_it_LU_decomposed_private;  
  unsigned int *row_permutation_private; 
};

template <typename SCALAR_TYPE>
SCALAR_TYPE & matrix<SCALAR_TYPE>::operator() (const unsigned int i,const unsigned int j) const
{
  return table[i][j];
}



template <typename SCALAR_TYPE>
const class matrix<SCALAR_TYPE>& matrix<SCALAR_TYPE>::operator = (const class matrix<SCALAR_TYPE> &A)
{ 
  is_it_LU_decomposed_private = A.is_it_LU_decomposed ();

  for (unsigned int i = 0 ; i < dim ; i++)
    row_permutation_private[i] = A.row_permutation(i);


  class matrix<SCALAR_TYPE> &M = *this;

  for (unsigned int i = 0 ; i < dim ; i++)
    for (unsigned int j = 0 ; j < dim ; j++)  
      M(i,j) = A(i,j);

  return M;
}




// Matrix conjugue
template <typename SCALAR_TYPE>
class matrix<SCALAR_TYPE> conj(const class matrix<SCALAR_TYPE> &M)
{
  class matrix<SCALAR_TYPE> result (M.dim);

  for (unsigned int i = 0 ; i < M.dim ; i++)
    for (unsigned int j = 0 ; j < M.dim ; j++)
      {result(j,i) = conj(M(i,j));}

  return result;
}



// "+" overloading.
template <typename SCALAR_TYPE>  
class matrix<SCALAR_TYPE> matrix<SCALAR_TYPE>::operator + (const class matrix<SCALAR_TYPE> &A) 
{ 
  const class matrix<SCALAR_TYPE> &M = *this;
  class matrix<SCALAR_TYPE> result(dim);

  for (unsigned int i = 0 ; i < dim ; i++)
    for (unsigned int j = 0 ; j < dim ; j++)  
      result(i,j) = M(i,j) + A(i,j); 

  return result;
}



// "+=" overloading.
template <typename SCALAR_TYPE>
class matrix<SCALAR_TYPE>& matrix<SCALAR_TYPE>::operator += (const class matrix<SCALAR_TYPE> &A)
{
  class matrix<SCALAR_TYPE> &M = *this;

  for (unsigned int i = 0 ; i < dim ; i++)
    for (unsigned int j = 0 ; j < dim ; j++)  
      M(i,j) += A(i,j); 

  return M;
}
  


  
// "-" overloading.  
template <typename SCALAR_TYPE>
class matrix<SCALAR_TYPE> matrix<SCALAR_TYPE>::operator - (const class matrix<SCALAR_TYPE> &A) 
{
  const class matrix<SCALAR_TYPE> &M = *this;
  class matrix<SCALAR_TYPE> result(dim);

  for (unsigned int i = 0 ; i < dim ; i++)
    for (unsigned int j = 0 ; j < dim ; j++)  
      result(i,j) = M(i,j) - A(i,j); 

  return result;
}

  
// "-=" overloading.
template <typename SCALAR_TYPE>
class matrix<SCALAR_TYPE>& matrix<SCALAR_TYPE>::operator -= (const class matrix<SCALAR_TYPE>& A)
{ 
  class matrix<SCALAR_TYPE> &M = *this;
  
  for (unsigned int i = 0 ; i < dim ; i++)
    for (unsigned int j = 0 ; j < dim ; j++)  
      M(i,j) -= A(i,j); 
  
  return M;
}
  


// "*" overloading : internal law.

template <typename SCALAR_TYPE> 
class matrix<SCALAR_TYPE> matrix<SCALAR_TYPE>::operator * (const class matrix<SCALAR_TYPE> &A) 
{
  const class matrix<SCALAR_TYPE> &M = *this;
  class matrix<SCALAR_TYPE> result(dim);

  for (unsigned int i = 0 ; i < dim ; i++)
    for (unsigned int j = 0 ; j < dim ; j++)
    {
      result(i,j) = 0;
      for (unsigned int k = 0 ; k < dim ; k++)
	result(i,j) += M(i,k)*A(k,j); 
    }

  return result;
}


// "*" overloading : external law.
// A*x is overloaded.
template <typename SCALAR_TYPE>  
class  matrix<SCALAR_TYPE> matrix<SCALAR_TYPE>::operator * (const SCALAR_TYPE x) const
{
  class matrix<SCALAR_TYPE> result(dim);
  const class matrix<SCALAR_TYPE> &M = *this;

  for (unsigned int i = 0 ; i < dim ; i++)
    for (unsigned int j = 0 ; j < dim ; j++)
      result(i,j) = x*M(i,j); 
    
  return result;
}



// "*" overloading : external law.  
// It is a "outside" function because one has to overload x*A, which cannot be done in the class.
template <typename SCALAR_TYPE>
class matrix<SCALAR_TYPE> operator * (const SCALAR_TYPE x,const class matrix<SCALAR_TYPE> &A) 
{
  class matrix<SCALAR_TYPE> result(A.dim,A.dim);
   
  for (unsigned int i = 0 ; i < A.dim ; i++)
    for (unsigned int j = 0 ; j < A.dim ; j++)
      result(i,j) = x*A(i,j); 
    
  return result;
}



  
// "*=" overloading : external law.
template <typename SCALAR_TYPE>
class matrix<SCALAR_TYPE>& matrix<SCALAR_TYPE>::operator *= (const SCALAR_TYPE x)
{
  class matrix<SCALAR_TYPE> &M = *this;

  for (unsigned int i = 0 ; i < dim ; i++)
    for (unsigned int j = 0 ; j < dim ; j++)
      M(i,j) *= x; 

 return M;
}




// "/" overloading : external law.
// A/x is overloaded.
template <typename SCALAR_TYPE>  
class matrix<SCALAR_TYPE> matrix<SCALAR_TYPE>::operator / (const SCALAR_TYPE x) 
{
  class matrix<SCALAR_TYPE> result(dim);  
  class matrix<SCALAR_TYPE> &M = *this;
 
  const SCALAR_TYPE one_over_x = 1.0/x;

  for (unsigned int i = 0 ; i < dim ; i++)
    for (unsigned int j = 0 ; j < dim ; j++)
      {
	result(i,j) = M(i,j)*one_over_x; 
      }

  return result;
}


  
// "/=" overloading : external law.
template <typename SCALAR_TYPE>
class matrix<SCALAR_TYPE>& matrix<SCALAR_TYPE>::operator /= (const SCALAR_TYPE x)
{
  const SCALAR_TYPE one_over_x = 1.0/x;
  class matrix<SCALAR_TYPE> &M = *this;

  for (unsigned int i = 0 ; i < dim ; i++)
    for (unsigned int j = 0 ; j < dim ; j++)
      M(i,j) *= one_over_x; 

 return M;
}

  


  


// Infinite norm.
template <typename SCALAR_TYPE>
double matrix<SCALAR_TYPE>::infinite_norm (void) const 
{
  const class matrix<SCALAR_TYPE> &M = *this;
  double try_norm,norm = 0.0;

  for (unsigned int i = 0 ; i < dim ; i++)
    for (unsigned int j = 0 ; j < dim ; j++)
      if (norm<(try_norm = abs (M(i,j))))
	norm = try_norm;
  
  return norm;
}
    
                                                                                                                                             
// trace.
template <typename SCALAR_TYPE>
complex <double> matrix<SCALAR_TYPE>::trace (void) const
{
  const class matrix<SCALAR_TYPE> &M = *this;
  complex <double> tr=0.0;
                                                                                                                                             
  for (unsigned int i = 0 ; i < dim ; i++)
	{
           tr+=M(i,i);
         }                                                                                                                                    
  return tr;
}
                                                                                                                                             
                                                                                                                                            
    
     
// Square matrix only.
template <typename SCALAR_TYPE>
void matrix<SCALAR_TYPE>::identity (void)
{
  if (dim !=dim) {cout<<"Identity is impossible"<<endl; exit(1);}

  class matrix<SCALAR_TYPE> &M = *this;
 
  for (unsigned int i = 0 ; i < dim ; i++)
    for (unsigned int j = 0 ; j < i ; j++)
       M(i,j) = M(j,i) = 0;
 
   for (unsigned int i = 0 ; i < dim ; i++) M(i,i) = 1;
}
 
   
template <typename SCALAR_TYPE>
void matrix<SCALAR_TYPE>::zero (void)
{
  class matrix<SCALAR_TYPE> &M = *this;
 
  for (unsigned int i = 0 ; i < dim ; i++)
    for (unsigned int j = 0 ; j < dim ; j++)
       M(i,j) = 0;
}




// Matrix transpose.
template <typename SCALAR_TYPE>
class matrix<SCALAR_TYPE> transpose_matrix (const class matrix<SCALAR_TYPE> &M)
{
  class matrix<SCALAR_TYPE> result (M.dim);

  for (unsigned int i = 0 ; i < M.dim ; i++)
    for (unsigned int j = 0 ; j < M.dim ; j++)
      result(j,i) = M(i,j);

  return result;
}




// Matrix to vectors for square matrices.
template <typename SCALAR_TYPE> 
void matrix<SCALAR_TYPE>::matrix_to_vectors (class vector_class<SCALAR_TYPE> **&eigenvectors)
{
  class matrix<SCALAR_TYPE> &M = *this;

  M = transpose_matrix (M);

  eigenvectors = new class vector_class<SCALAR_TYPE> * [dim];

  for (unsigned int i = 0 ; i < dim ; i++)
  {
    eigenvectors[i] = new class vector_class<SCALAR_TYPE> (dim); 
    class vector_class<SCALAR_TYPE> &Vi = *(eigenvectors[i]);

    for (unsigned int j = 0 ; j < dim ; j++)
      Vi(j) = M(i,j);

    delete [] table[i];
  }

  delete [] table;
  deleted = true;
}



// Swapping of two row vectors of the matrix.
// ------------------------------------------

template <typename SCALAR_TYPE>
void matrix<SCALAR_TYPE>::swap_vectors (const unsigned int i,const unsigned int j)
{
  if (i == j) return;

  class vector_class<SCALAR_TYPE> *const temp_eig_vec = table[i];
  table[i] = table[j];
  table[j] = temp_eig_vec;
}


// V=AX.

template <typename SCALAR_TYPE>
class vector_class<SCALAR_TYPE> operator * (const class matrix<SCALAR_TYPE> &A,
					    const class vector_class<SCALAR_TYPE> &X)
{
  class vector_class<SCALAR_TYPE> V(A.dim);  
 
  for (unsigned int i = 0 ; i < A.dim ; i++)
  {
    V(i) = 0.0;

    for (unsigned int j = 0 ; j < A.dim ; j++)
      {V(i) += A(i,j)*X(j); }
  }

  return V;
}




// tensor product of two vectors => matrix.

template <typename SCALAR_TYPE> 
class matrix<SCALAR_TYPE> tensor_product (const class vector_class<SCALAR_TYPE> &A,
					  const class vector_class<SCALAR_TYPE> &B) 
{
  class matrix<SCALAR_TYPE> result(A.dim);

  for (unsigned int i = 0 ; i < A.dim ; i++)
    for (unsigned int j = 0 ; j < A.dim ; j++)
      result(i,j) = A(i)*B(j);

  return result;
}




template <typename SCALAR_TYPE>
 ostream & operator << ( ostream &os,
			     const class matrix<SCALAR_TYPE> &M)
{
  for (unsigned int i = 0 ; i < M.dim ; i++)
  {
    for (unsigned int j = 0 ; j < M.dim ; j++)
      os<<M(i,j)<<" ";

    os<<endl;
  }
  return os;
}


template <typename SCALAR_TYPE>
 istream & operator >> ( istream &is,
			     class matrix<SCALAR_TYPE> &M)
{
  for (unsigned int i = 0 ; i < M.dim ; i++)
    for (unsigned int j = 0 ; j < M.dim ; j++)
      is>>M(i,j);

  return is;
}



// LU decomposition of the matrix.
// -------------------------------
// It is used to solve the linear system AX = B.

template <typename SCALAR_TYPE>
void matrix <SCALAR_TYPE>::LU_decompose (void)
{
  class matrix <SCALAR_TYPE> &A = *this;
  class vector_class<double> V(dim);
  unsigned int imax;
  double temp1;
  for (unsigned int i = 0 ; i < dim ; i++)
  {
    double big = 0.0;
    for (unsigned int j = 0 ; j < dim ; j++)
      if ((temp1 = inf_norm (A(i,j))) > big)
	big = temp1;

    V(i) = 1.0/big;
  }

  for (unsigned int j = 0 ; j < dim ; j++)
  {
    for (unsigned int i = 0 ; i < j ; i++)
    {
      SCALAR_TYPE sum = A(i,j);
      for (unsigned int k = 0 ; k < i ; k++) sum -= A(i,k)*A(k,j);

      A(i,j) = sum;
    }

    double big = 0.0;
    for (unsigned int i = j ; i < dim ; i++)
    {
      SCALAR_TYPE sum = A(i,j);
      for (unsigned int k = 0 ; k < j ; k++) sum -= A(i,k)*A(k,j);

      A(i,j) = sum;

      if ((temp1 = V(i)*inf_norm (sum)) >= big)
      {
	big = temp1;
	imax = i;
      }
    }
    
    if (j != imax)
    {
      for (unsigned int k = 0 ; k < dim ; k++)
      {
	const SCALAR_TYPE temp2 = A(imax,k);
	A(imax,k) = A(j,k);
	A(j,k) = temp2;
      }
  
      V(imax) = V(j);
    }
    
    row_permutation(j) = imax;
  
    if (A(j,j) == 0.0) A(j,j) = precision;
    
    if (j != dim-1)
    {
      const SCALAR_TYPE temp2 = 1.0/A(j,j);

      for (unsigned int i = j+1 ; i < dim ;i++)
	A(i,j) *= temp2;
    }
  }

  is_it_LU_decomposed_private = true;
}
