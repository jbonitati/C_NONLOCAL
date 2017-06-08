// Namespace of routines diagonalizing a complex symmetric matrix.
// ---------------------------------------------------------------
// One can have tridiagonal or general matrices.
// One can have eigenvalues (semi diagonalization) or eigenvalues and eigenvectors (diagonalization).
// The method is to use the Householder method to put the matrix under tridiagonal form, (if not already)
// and after use the QL diagonalization procedure with the tridiagonal matrix.
// See Numerical Recipes for more details.
// Eigenvectors are normalized with the Berggren norm.
// Eigenvectors and eigenvalues are sorted according to the real part of the eigenvalues.

namespace Householder_QL
{
  template <typename SCALAR_TYPE> 
  void all_eigenpairs (class matrix<SCALAR_TYPE> &A,
		       class array<SCALAR_TYPE> &eigenvalues);

  template <typename SCALAR_TYPE> 
  void all_eigenpairs_inverse_iteration (class matrix<SCALAR_TYPE> &A,
					 class array<SCALAR_TYPE> &eigenvalues,
					 const double prec);

  template <typename SCALAR_TYPE> 
  void all_eigenpairs (const class array<SCALAR_TYPE> &diag,
		       const class array<SCALAR_TYPE> &off_diag,
		       class matrix<SCALAR_TYPE> &A,
		       class array<SCALAR_TYPE> &eigenvalues);

  template <typename SCALAR_TYPE> 
  void all_eigenvalues (class matrix<SCALAR_TYPE> &A,
			class array<SCALAR_TYPE> &eigenvalues);

  template <typename SCALAR_TYPE> 
  void all_eigenvalues (const class array<SCALAR_TYPE> &diag,
			const class array<SCALAR_TYPE> &off_diag,
			class array<SCALAR_TYPE> &eigenvalues);

  template <typename SCALAR_TYPE> 
  void eigen_sort (const int low,const int high,
		   class array<SCALAR_TYPE> &eigenvalues,
		   class matrix<SCALAR_TYPE> *const A_ptr);
  
  template <typename SCALAR_TYPE> 
  void tridiag (const bool is_it_only_eigenvalues,
		class matrix<SCALAR_TYPE> &A,
		class array<SCALAR_TYPE> &diag,
		class array<SCALAR_TYPE> &off_diag);

  template <typename SCALAR_TYPE> 
  void QL_diag (class matrix<SCALAR_TYPE> *const A_ptr,
		class array<SCALAR_TYPE> &diag,
		class array<SCALAR_TYPE> &off_diag);
}







// Sort of eigenvalues and/or eigenvectors.
// ----------------------------------------
// Quicksort is used to sort eigenvalues and/or eigenvectors.
// The sort is made with the real part of eigenvalues.
//
// Variables :
// -----------
// low,high : boundaries used in Quicksort.
// pivot_eigenvalue : Quicksort pivot = real (eigenvalues((low + high)/2))
// eigenvalues : set of eigenvalues to sort.
// A_ptr : pointer on the matrix of eigenvectors. 
//         It is 0 if one wants only eigenvalues.

template <typename SCALAR_TYPE> 
void Householder_QL::eigen_sort (const int low,const int high,
				 class array<SCALAR_TYPE> &eigenvalues,
				 class matrix<SCALAR_TYPE> *const A_ptr)
{
  int i_sort = low;
  int j_sort = high;
  const double pivot_eigenvalue = real (eigenvalues((i_sort + j_sort)/2));

  do
  {
    while ((real (eigenvalues(i_sort)) < pivot_eigenvalue) && (abs (real (eigenvalues(i_sort)) - pivot_eigenvalue) > precision))
      i_sort++;

    while ((real (eigenvalues(j_sort)) > pivot_eigenvalue) && (abs (real (eigenvalues(j_sort)) - pivot_eigenvalue) > precision))
      j_sort--;
     
    if (i_sort <= j_sort)
    {
      const SCALAR_TYPE temp_eigenvalue = eigenvalues(i_sort);
      eigenvalues(i_sort) = eigenvalues(j_sort);
      eigenvalues(j_sort) = temp_eigenvalue;

      if (A_ptr != 0) A_ptr->swap_vectors (i_sort,j_sort);
     
      i_sort++;
      j_sort--;
    }
  }
  while (i_sort <= j_sort);
  
  if (low < j_sort) eigen_sort (low,j_sort,eigenvalues,A_ptr);

  if (i_sort < high) eigen_sort (i_sort,high,eigenvalues,A_ptr); 
}






// Tridiagonalization of a general summetric matrix.
// -------------------------------------------------
// The Householder method is used here.
// It is generalized to the complex case.
// The routine is an adaptation of the freely available FORTRAN Netlib routine tred2. 
//
// Variables :
// -----------
// A : matrix to diagonalize. 
//     It contains at the end the matrix of passage from the total matrix to the tridiagonal matrix.
// diag : It will contain at the end the diagonal matrix elements of the tridiagonal matrix.
// off_diag : It will contain at the end the off diagonal matrix elements of the tridiagonal matrix.


template <typename SCALAR_TYPE> 
void Householder_QL::tridiag (const bool is_it_only_eigenvalues,
			      class matrix<SCALAR_TYPE> &A,
			      class array<SCALAR_TYPE> &diag,
			      class array<SCALAR_TYPE> &off_diag)
{
  const unsigned int N = A.dim;
  for (unsigned int i = 0 ; i < N ; i++) diag(i) = A(i,N-1);

  for (unsigned int i = N-1 ; i >= 1 ; i--)
  {
    SCALAR_TYPE h = 0.0;

    if (i == 1) { off_diag(i) = diag(i-1); for (unsigned int j = 0 ; j < i ; j++) { diag(j) = A(j,i-1); A(i,j) = A(j,i) = 0.0; }} 
    else
    {
      double scale = 0.0;
      for (unsigned int k = 0 ; k < i ; k++) scale += inf_norm (diag(k));

      if (scale == 0.0) { off_diag(i) = diag(i-1); for (unsigned int j = 0 ; j < i ; j++) { diag(j) = A(j,i-1); A(i,j) = A(j,i) = 0.0; }}
      else
      {
	for (unsigned int k = 0 ; k < i ; k++) { diag(k) /= scale; h += diag(k)*diag(k); }
	const SCALAR_TYPE sqrt_h = sqrt (h),f = diag(i-1),g = (inf_norm (f - sqrt_h) > inf_norm (f + sqrt_h)) ? (sqrt_h) : (-sqrt_h);
	off_diag(i) = g*scale; h -= f*g; diag(i-1) = f-g;
	
	for (unsigned int j = 0 ; j < i ; j++) off_diag(j) = 0.0;

	for (unsigned int j = 0 ; j < i ; j++) 
	{
	  A(i,j) = diag(j);
	  const SCALAR_TYPE ff = diag(j);
	  SCALAR_TYPE gg = off_diag(j) + A(j,j)*ff;
	  for (unsigned int k = j+1 ; k < i ; k++) { gg += A(j,k)*diag(k); off_diag(k) += A(j,k)*ff; }
	  off_diag(j) = gg; 
	}

	SCALAR_TYPE ff = 0.0;
	for (unsigned int j = 0 ; j < i ; j++) { off_diag(j) /= h; ff += off_diag(j)*diag(j); }
	const SCALAR_TYPE hh = ff/(h + h);
	for (unsigned int j = 0 ; j < i ; j++) off_diag(j) -= hh*diag(j);

	for (unsigned int j = 0 ; j < i ; j++) 
	{
	  const SCALAR_TYPE fff = diag(j),ggg = off_diag(j);
	  for (unsigned int k = j ; k < i ; k++) A(j,k) -= fff*off_diag(k) + ggg*diag(k);
	  diag(j) = A(j,i-1);
	  A(j,i) = 0.0;
	}
      }
    }
    diag(i) = h;
  }

  off_diag(0) = 0.0;

  if (is_it_only_eigenvalues) 
    for (unsigned int i = 0 ; i < N ; i++) diag(i) = A(i,i);
  else
  {
    for (unsigned int i = 1 ; i < N ; i++)
    {
      const SCALAR_TYPE h = diag(i);
      A(i-1,N-1) = A(i-1,i-1); A(i-1,i-1) = 1.0;

      if (h != 0.0)
      {
	for (unsigned int k = 0 ; k < i ; k++) diag(k) = A(i,k)/h;

	for (unsigned int j = 0 ; j < i ; j++)
	{
	  SCALAR_TYPE g = 0.0;
	  for (unsigned int k = 0 ; k < i ; k++) g += A(i,k)*A(j,k);
	  for (unsigned int k = 0 ; k < i ; k++) A(j,k) -= g*diag(k);
	}
      }

      for (unsigned int k = 0 ; k < i ; k++) A(i,k) = 0.0;
    }

    for (unsigned int i = 0 ; i < N ; i++) { diag(i) = A(i,N-1); A(i,N-1) = (i == N-1) ? (1.0) : (0.0); }
  }
}








// Diagonalization of a tridiagonal symmetric matrix.
// --------------------------------------------------
// The QL diagonalization with implicit shifts method is used here.
// This routine is an adaptation of the freely available FORTRAN Netlib routine imtql2.
// It is generalized to the complex case.
// See Numerical Recipes for more details about the method.
//
// Variables :
// -----------
// A_ptr : Pointer to a class matrix.
//         This pointer is sero if one does not want the eigenvectors.
//         If not, it is the matrix of passage from the total matrix to diagonalize to the tridiagonal matrix.
//         It will contain at the end the normalized eigenvectors of the total matrix. 
// A,dummy : reference on the A matrix if A_ptr != 0, or reference on dummy if A_ptr = 0. dummy is a dummy matrix with one element. 
// diag : Array of the diagonal matrix elements of the tridiagonal matrix.
//        It will contain at the end the eigenvalues of the tridiagonal matrix.
// off_diag : Array of the off diagonal matrix elements of the tridiagonal matrix.
//            Its values are meaningless at the end.

template <typename SCALAR_TYPE> 
void Householder_QL::QL_diag (class matrix<SCALAR_TYPE> *const A_ptr,
			      class array<SCALAR_TYPE> &diag,
			      class array<SCALAR_TYPE> &off_diag)
{   
  const bool is_it_only_eigenvalues = (A_ptr == 0);
  const class matrix<SCALAR_TYPE> dummy(1),&A = (!is_it_only_eigenvalues) ? (*A_ptr) : (dummy);
  const unsigned int N = diag.dim(0);

  for (unsigned int l = 0 ; l < N ; l++)
  { 
    unsigned int count = 0, m = N;

    while (m != l)
    {
      m = l;
      while ((m < N-1) && (inf_norm (off_diag(m))/(inf_norm (diag(m)) + inf_norm (diag(m+1))) >= 1E-13)) m++;
      
      if (m != l)
      {
	if (count++ == 30) {cout<<"No convergence in QL diagonalization."<<endl; exit (1);}
	if ((!isfinite (diag(m))) || (!isfinite (off_diag(m)))) {cout<<"Infinite numbers encountered during QL diagonalization."<<endl; exit (1);}
      
	const SCALAR_TYPE e = (diag(l+1) - diag(l))/(2.0*off_diag(l)),he = hypot (1.0,e),e_p_he = e + he,e_m_he = e - he;
	SCALAR_TYPE g = (inf_norm (e_p_he) > inf_norm (e_m_he)) ? (diag(m) - diag(l) + off_diag(l)/e_p_he) : (diag(m) - diag(l) + off_diag(l)/e_m_he);
	
	unsigned int i = m-1;
	SCALAR_TYPE r = 1.0,s = 1.0,c = 1.0,p = 0.0;
	while ((i >= l) && (i < N) && (r != 0.0))
	{
	  const SCALAR_TYPE f = s*off_diag(i),b = c*off_diag(i);
	  off_diag(i+1) = r = hypot (f,g);
	  
	  if (r != 0.0)
	  {
	    s = f/r; c = g/r; g = diag(i+1) - p; r = (diag(i) - g)*s + 2.0*c*b; p = s*r; diag(i+1) = g + p; g = c*r - b;
	    if (!is_it_only_eigenvalues) 
	      for (unsigned int k = 0 ; k < N ; k++) { const SCALAR_TYPE d = A(i+1,k); A(i+1,k) = s*A(i,k) + c*d; A(i,k) = c*A(i,k) - s*d; }
	    i--;
	  }
	  else { diag(i+1) -= p; off_diag(m) = 0.0; }
	}
	if ((r != 0.0) || (i < l)) { diag(l) -= p; off_diag(l) = g; off_diag(m) = 0.0;}       
      }}}
}




// Diagonalization of a general matrix.
// ------------------------------------
// Sorted eigenvalues and eigenvectors with the Householder + QL method.
// The sort is made with the real part of eigenvalues.
//
// Variables :
// -----------
// A : Matrix to diagonalize.
//     It will contain at the end the sorted normalized eigenvectors of the total matrix.
// eigenvalues : It will contain at the end the sorted eigenvalues of the tridiagonal matrix.
// off_diag : Array of the off diagonal matrix elements of the tridiagonal matrix.

template <typename SCALAR_TYPE> 
void Householder_QL::all_eigenpairs (class matrix<SCALAR_TYPE> &A,
				     class array<SCALAR_TYPE> &eigenvalues)
{ 
  const unsigned int N = A.dim;
  
  if (N == 1) 
  {
    eigenvalues(0) = A(0,0);
    A(0,0) = 1.0;
    return;
  }

  class array<SCALAR_TYPE> off_diag(N);

  tridiag (false,A,eigenvalues,off_diag);

  for (unsigned int i = 1 ; i < N ; i++) off_diag(i-1) = off_diag(i);
  off_diag(N-1) = 0.0;

  QL_diag (&A,eigenvalues,off_diag);
  eigen_sort (0,N-1,eigenvalues,&A);
}






// Diagonalization of a tridiagonal matrix.
// ----------------------------------------
// Sorted eigenvalues and eigenvectors with the QL method.
// The sort is made with the real part of eigenvalues.
//
// Variables :
// -----------
// diag : Array of the diagonal matrix elements of the tridiagonal matrix.
// off_diag : Array of the off diagonal matrix elements of the tridiagonal matrix.
// off_diag_temp : Copy of off_diag. It is used as its values change and become meaningless at the end of the calculation.
// A : It will contain at the end the sorted normalized eigenvectors of the tridiagonal matrix.
// eigenvalues : It will contain at the end the sorted eigenvalues of the tridiagonal matrix.

template <typename SCALAR_TYPE> 
void Householder_QL::all_eigenpairs (const class array<SCALAR_TYPE> &diag,
				     const class array<SCALAR_TYPE> &off_diag,
				     class matrix<SCALAR_TYPE> &A,
				     class array<SCALAR_TYPE> &eigenvalues)
{
  const unsigned int N = A.dim;

  if (N == 1) 
  {
    eigenvalues(0) = diag(0);
    A(0,0) = 1.0;
    return;
  }

  class array<SCALAR_TYPE> off_diag_temp (N);

  for (unsigned int i = 0 ; i < N ; i++) off_diag_temp(i) = off_diag(i);
  off_diag_temp(N-1) = 0.0;

  for (unsigned int i = 0 ; i < N ; i++) eigenvalues(i) = diag(i);
 
  A.identity ();
  QL_diag (&A,eigenvalues,off_diag_temp);

  eigen_sort (0,N-1,eigenvalues,&A);
}


template <typename SCALAR_TYPE> 
void Householder_QL::all_eigenpairs_inverse_iteration (class matrix<SCALAR_TYPE> &A,
						       class array<SCALAR_TYPE> &eigenvalues,
						       const double prec)
{ 
  const double A_inf = A.infinite_norm ();
  A /= A_inf;

  const class matrix<SCALAR_TYPE> AA = A;

  const unsigned int N = A.dim;

  for (unsigned int i = 0 ; i < N ; i++) 
    for (unsigned int j = 0 ; j < i  ; j++) 
      if (inf_norm (A(i,j)) < prec) A(i,j) = A(j,i) = 0.0;

  all_eigenpairs (A,eigenvalues);

  for (unsigned int i = 0 ; i < N ; i++)
    for (unsigned int loop = 0 ; loop < 2 ; loop++)
    {
      class matrix<SCALAR_TYPE> A_eig_i = AA;
      for (unsigned int ii = 0 ; ii < N ; ii++) A_eig_i(ii,ii) -= eigenvalues(i);

      A.eigenvector(i) = linear_system_solution (A_eig_i,A.eigenvector(i));
      A.eigenvector(i).normalization ();
      eigenvalues(i) = A.eigenvector(i)*(AA*A.eigenvector(i));
    }

  for (unsigned int i = 0 ; i < N ; i++) eigenvalues(i) *= A_inf;
}



// Semi diagonalization of a general matrix.
// -----------------------------------------
// Sorted eigenvalues with the Householder + QL method.
// The sort is made with the real part of eigenvalues.
//
// Variables :
// -----------
// A : Matrix to diagonalize. It is erased at the end of the calculation.
// eigenvalues : It will contain at the end the sorted eigenvalues of the tridiagonal matrix.
// off_diag : Array of the off diagonal matrix elements of the tridiagonal matrix.

template <typename SCALAR_TYPE> 
void Householder_QL::all_eigenvalues (class matrix<SCALAR_TYPE> &A,
				      class array<SCALAR_TYPE> &eigenvalues)
{
  const unsigned int N = A.dim;

  if (N == 1) 
  {
    eigenvalues(0) = A(0,0);
    return;
  }

  class array<SCALAR_TYPE> off_diag(N);

  tridiag (true,A,eigenvalues,off_diag);

  for (unsigned int i = 1 ; i < N ; i++) off_diag(i-1) = off_diag(i);
  off_diag(N-1) = 0.0;

  class matrix<SCALAR_TYPE> *const NO_MATRIX = 0;

  QL_diag (NO_MATRIX,eigenvalues,off_diag);

  eigen_sort (0,N-1,eigenvalues,NO_MATRIX);
}






// Semi diagonalization of a tridiagonal matrix.
// ---------------------------------------------
// Sorted eigenvalues with the QL method.
// The sort is made with the real part of eigenvalues.
//
// Variables :
// -----------
// diag : Array of the diagonal matrix elements of the tridiagonal matrix.
// off_diag : Array of the off diagonal matrix elements of the tridiagonal matrix.
// off_diag_temp : Copy of off_diag. It is used as its values change and become meaningless at the end of the calculation.
// eigenvalues : It will contain at the end the sorted eigenvalues of the tridiagonal matrix.

template <typename SCALAR_TYPE> 
void Householder_QL::all_eigenvalues (const class array<SCALAR_TYPE> &diag,
				      const class array<SCALAR_TYPE> &off_diag,
				      class array<SCALAR_TYPE> &eigenvalues)
{
  const unsigned int N = diag.dim (0);

  if (N == 1) 
  {
    eigenvalues(0) = diag(0,0);
    return;
  }

  class array<SCALAR_TYPE> off_diag_temp (N);

  for (unsigned int i = 0 ; i < N-1 ; i++) off_diag_temp(i) = off_diag(i);
  off_diag_temp(N-1) = 0.0;

  class matrix<SCALAR_TYPE> *const NO_MATRIX = 0;

  for (unsigned int i = 0 ; i < N ; i++) eigenvalues(i) = diag(i);

  QL_diag (NO_MATRIX,eigenvalues,off_diag_temp);

  eigen_sort (0,N-1,eigenvalues,NO_MATRIX);
}
