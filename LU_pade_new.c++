// Linear system solver.
// ---------------------
// When the matrix A has been LU decomposed, one uses this routine to get X in AX=B.

template <typename SCALAR_TYPE>
class vector_class<SCALAR_TYPE> linear_system_solution (class matrix<SCALAR_TYPE> &A,
							const class vector_class<SCALAR_TYPE> &B)
{
  if (A.is_it_LU_decomposed () == false) 
    A.LU_decompose ();

  const unsigned int N = A.dim; 
  unsigned int ii = N;

  SCALAR_TYPE sum;
  class vector_class<SCALAR_TYPE> V = B;

  for (unsigned int i = 0 ; i < N ; i++)
  {
    const unsigned int ip = A.row_permutation(i);
    sum = V(ip);
    V(ip) = V(i);

    if (ii != N)
      for (int j = ii ; j <= i-1 ; j++)
	sum -= A(i,j)*V(j);
    else if (sum != 0.0)
      ii = i;

    V(i) = sum;
  }

  for (unsigned int i = N-1 ; i <= N-1 ; i--)
  {
    sum = V(i);

    for (unsigned int j = i+1 ; j < N ; j++)
      sum -= A(i,j)*V(j);
    
    V(i) = sum/A(i,i);
  }

  return V;
} 



