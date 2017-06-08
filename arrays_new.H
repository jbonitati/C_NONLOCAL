template <class C>
class array_dim1
{
public:
  array_dim1 (const unsigned int d0_c) 
    : d0(d0_c)
  {
    if (d0 == 0)
      {
	cout<<"Dimension is zero."<<endl;
	abort();
      }
    
    table = new C [d0];
  }
  
  ~array_dim1 (void)
  {
    delete [] table;
  }
  
  C *table;
  
  C & operator () (const unsigned int i0) const;
  
  const unsigned int d0;
};

template <class C>
C & array_dim1<C>::operator () (const unsigned int i0) const
{
  if (i0 >= d0)
    {
      cout<<"Problem with i0="<<i0<<" >= d0="<<d0<<"."<<endl;
      abort ();
    }
  else
    return table[i0];
}



template <class C>
class array_dim2
{
public:
  
  array_dim2 (const unsigned int d0_c,const unsigned int d1_c) 
    : d0(d0_c), d1(d1_c)
  {
    if ( (d0 == 0) || (d1 == 0))
      {
	cout<<"Dimension is zero."<<endl;
	abort();
      }
    
    table = new C [d0*d1];
  }

  ~array_dim2 (void)
  {
    delete [] table;
  }
  
  C *table;
  
  C & operator () (const unsigned int i0,const unsigned int i1) const;
  
  const unsigned int d0,d1;
};


template <class C>
C & array_dim2<C>::operator () (const unsigned int i0,const unsigned int i1) const
{

  if (d1 == 0)
    {
      cout<<"1 component too much."<<endl;
      abort ();
    }
  
  else if (i0 >= d0)
    {
      cout<<"Problem with i0="<<i0<<" >= d0="<<d0<<"."<<endl;
      abort ();
    }
  
  else if (i1 >= d1)
    {
      cout<<"Problem with i1="<<i1<<" >= d1="<<d1<<"."<<endl;
      abort ();
    }
  
  else
    return table[i1 + d1*i0];
}
