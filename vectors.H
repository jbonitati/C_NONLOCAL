template <typename SCALAR_TYPE>
class vector_class
{
public :
  vector_class (const unsigned int dim_CREATE) : dim (dim_CREATE)
    {
      list = new SCALAR_TYPE [dim];

      for (unsigned int i = 0 ; i < dim ; i++)
	list[i] = 0;
    }
  
  vector_class (const class vector_class &V) : dim (V.dim)
    { 
      list = new SCALAR_TYPE [dim];
      
      for (unsigned int i = 0 ; i < dim ; i++)
	list[i] = V.list[i];  
    }

 ~vector_class (void)
    {
      delete [] list;
    }

  class vector_class& operator += (const class vector_class &);
  class vector_class& operator -= (const class vector_class &);
  class vector_class& operator *= (const SCALAR_TYPE );
  class vector_class& operator /= (const SCALAR_TYPE );

  const class vector_class& operator = (const vector_class &);
  class vector_class operator + (const class vector_class &) const;
  class vector_class operator - (const class vector_class &) const;
  class vector_class operator * (const SCALAR_TYPE ) const; 
  class vector_class operator / (const SCALAR_TYPE ) const; 

  SCALAR_TYPE & operator() (const unsigned int) const;

  double infinite_norm (void) const;
  complex<double> coef_max(void) const;
  int index_max(void) const ;


  SCALAR_TYPE dot_product_norm (void) const;
  void normalization (void);
  void zero (void); 
  void random_vector (void);

  const unsigned int dim;  

private:
  SCALAR_TYPE *list; 
};

template <typename SCALAR_TYPE>
SCALAR_TYPE & vector_class<SCALAR_TYPE>::operator() (const unsigned int i) const
{
  return list[i];
}




// The matrices and vector_classs must have the same dim. No test is made.


template <typename SCALAR_TYPE>
const class vector_class<SCALAR_TYPE>& vector_class<SCALAR_TYPE>::operator = (const class vector_class<SCALAR_TYPE> &V)
{
  class vector_class<SCALAR_TYPE> &W = *this;

  for (unsigned int i = 0 ; i < dim ; i++)
    W(i) = V(i); 

  return W;
}

// "+" overloading.
template <typename SCALAR_TYPE>  
class vector_class<SCALAR_TYPE> vector_class<SCALAR_TYPE>::operator + (const class vector_class<SCALAR_TYPE> &V) const
{
  const class vector_class<SCALAR_TYPE> &W = *this;
  class vector_class<SCALAR_TYPE> result(dim);

  for (unsigned int i = 0 ; i < dim ; i++)
    result(i) = W(i) + V(i); 

  return result;
}

// "+=" overloading.
template <typename SCALAR_TYPE>
class vector_class<SCALAR_TYPE>& vector_class<SCALAR_TYPE>::operator += (const class vector_class<SCALAR_TYPE> &V)
{
  class vector_class<SCALAR_TYPE> &W = *this;

  for (unsigned int i = 0 ; i < dim ; i++)
    W(i) += V(i); 

  return W;
}
  


  
// "-" overloading.  
template <typename SCALAR_TYPE>
class vector_class<SCALAR_TYPE> vector_class<SCALAR_TYPE>::operator - (const class vector_class<SCALAR_TYPE> &V) const
{
  const class vector_class<SCALAR_TYPE> &W = *this;
  class vector_class<SCALAR_TYPE> result(dim);

  for (unsigned int i = 0 ; i < dim ; i++)
    result(i) = W(i) - V(i); 

  return result;
}

  
// "-=" overloading.
template <typename SCALAR_TYPE>
class vector_class<SCALAR_TYPE>& vector_class<SCALAR_TYPE>::operator -= (const class vector_class<SCALAR_TYPE>& V)
{
  class vector_class<SCALAR_TYPE> &W = *this;

  for (unsigned int i=0;i<dim;i++)
    W(i) -= V(i);
  
 return W;
}
  

  

// "*" overloading : external law.
// V*x is overloaded.
template <typename SCALAR_TYPE>  
class vector_class<SCALAR_TYPE> vector_class<SCALAR_TYPE>::operator * (const SCALAR_TYPE x) const
{
  const class vector_class<SCALAR_TYPE> &W = *this;
  class vector_class<SCALAR_TYPE> result(dim);

  for (unsigned int i=0;i<dim;i++)
      result(i) = x*W(i); 
    
  return result;
}

// "*" overloading : external law.  
// It is a "outside" function because one has to overload x*V, which cannot be done in the class.
template <typename SCALAR_TYPE>
class vector_class<SCALAR_TYPE> operator * (const SCALAR_TYPE x,const class vector_class<SCALAR_TYPE> &V) 
{
  class vector_class<SCALAR_TYPE> result(V.dim);

  for (unsigned int i = 0 ; i < V.dim ; i++)
    result(i) = x*V(i); 
    
  return result;
}
  
// "*=" overloading : (external law)
template <typename SCALAR_TYPE>
class vector_class<SCALAR_TYPE>& vector_class<SCALAR_TYPE>::operator *= (const SCALAR_TYPE x)
{
  class vector_class<SCALAR_TYPE> &W = *this;

  for (unsigned int i = 0 ; i < dim ; i++)
    W(i) *= x; 

 return W;
}

// "/" overloading : external law.
// V/x is overloaded.
template <typename SCALAR_TYPE>  
class vector_class<SCALAR_TYPE> vector_class<SCALAR_TYPE>::operator / (const SCALAR_TYPE x) const
{
  const class vector_class<SCALAR_TYPE> &W = *this;
  class vector_class<SCALAR_TYPE> result(dim);
  SCALAR_TYPE one_over_x = 1.0/x;

  for (unsigned int i = 0 ; i < dim ; i++)
    result(i) = W(i)*one_over_x; 
    
  return result;
}

  
// "/=" overloading : (external law)
template <typename SCALAR_TYPE>
class vector_class<SCALAR_TYPE>& vector_class<SCALAR_TYPE>::operator /= (const SCALAR_TYPE x)
{
  class vector_class<SCALAR_TYPE> &W = *this;
  SCALAR_TYPE one_over_x = 1.0/x;	

  for (unsigned int i = 0 ; i < dim ; i++)
    W(i) *= one_over_x; 

 return W;
}


// Infinite norm.
template <typename SCALAR_TYPE> double vector_class<SCALAR_TYPE>::infinite_norm (void) const
{
  const class vector_class<SCALAR_TYPE> &W = *this;
  double try_norm,norm = 0.0;

  for (unsigned int i = 0 ; i < dim ; i++)
    if (norm < (try_norm = abs (W(i))))
      norm = try_norm;
  
  return norm;
}

//coef max
template <typename SCALAR_TYPE> complex<double> vector_class<SCALAR_TYPE>::coef_max (void) const
{
  const class vector_class<SCALAR_TYPE> &W = *this;
  complex<double> try_z;
  complex<double> coef = 0.0;

  for (unsigned int i = 0 ; i < dim ; i++)
    { 
      try_z =W(i);
      
      if ( abs(coef) < abs(try_z)     )
	{
	  coef = try_z;
	}
    }
  
  return coef;
}


//index max
template <typename SCALAR_TYPE> int vector_class<SCALAR_TYPE>::index_max (void) const
{
  const class vector_class<SCALAR_TYPE> &W = *this;
  complex<double> try_z;
  complex<double> coef = 0.0;
  int index=0;

  for (unsigned int i = 0 ; i < dim ; i++)
    { 
      try_z =W(i);
      
      if ( abs(coef) < abs(try_z)     )
	{
	  coef = try_z;
	  index=i;
	}
    }
  
  return index;
}







//  norm 2 : not the hermitian one for complex, but the Berggren one.
template <typename SCALAR_TYPE>
SCALAR_TYPE vector_class<SCALAR_TYPE>::dot_product_norm (void) const
{
  const class vector_class<SCALAR_TYPE> &W = *this;
  SCALAR_TYPE norm_square = 0.0;

  for (unsigned int i = 0 ; i < dim ; i++)
    norm_square += W(i)*W(i);
  
  return sqrt (norm_square);
}


// Dot product : not the hermitian one for complex, but the Berggren one.
template <typename SCALAR_TYPE>
SCALAR_TYPE operator * (const vector_class<SCALAR_TYPE> &X,const class vector_class<SCALAR_TYPE> &Y)
{
  SCALAR_TYPE product = 0.0;

  for (unsigned int i = 0 ; i < X.dim ; i++)
    product += X(i)*Y(i);
  
  return product;
}

//  normalization with the previous norm.
template <typename SCALAR_TYPE>
void vector_class<SCALAR_TYPE>::normalization (void)
{
  class vector_class<SCALAR_TYPE> &W = *this;

  W /= dot_product_norm ();
}  


//  V=0.
template <typename SCALAR_TYPE>
void vector_class<SCALAR_TYPE>::zero (void)
{
  class vector_class<SCALAR_TYPE> &W = *this;
  
  for (unsigned int i = 0 ; i < dim ; i++)  
    W(i) = 0.0;
}  

//  V=random.
template <typename SCALAR_TYPE>
void vector_class<SCALAR_TYPE>::random_vector (void)
{
  class vector_class<SCALAR_TYPE> &W = *this;

  for (unsigned int i = 0 ; i < dim ; i++)  
    W(i) = drand48 ();
}

  
template <typename SCALAR_TYPE>
 ostream & operator << ( ostream &os,
			     const class vector_class<SCALAR_TYPE> &V)
{
  for (unsigned int i = 0 ; i < V.dim ; i++)
    os<<V(i)<<" ";
  
  return os;
}


template <typename SCALAR_TYPE>
 istream & operator >> ( istream &is,
			     class vector_class<SCALAR_TYPE> &V)
{
  for (unsigned int i = 0 ; i < V.dim ; i++)
    is>>V(i);

  return is;
}
