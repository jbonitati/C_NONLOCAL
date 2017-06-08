
template <class C>
class array_large
{
public:

  array_large (const uint64_t d0_c) 
    : d0(d0_c)
  {
    if (d0 == 0)
    {
      cout<<"Dimension is zero."<<endl;
      abort();
    }
    

    table = new C [d0];
  }





  ~array_large (void)
  {
    delete [] table;
  }
  
  C *table;

 
  C & operator [] (const uint64_t i0) const;

  C & operator () (const uint64_t i0) const;
  


  const uint64_t d0;
};




template <class C>
C & array_large<C>::operator [] (const  uint64_t i0) const
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
C & array_large<C>::operator () (const uint64_t i0) const
{
  
  if (i0 >= d0)
  {
    cout<<"Problem with i0="<<i0<<" >= d0="<<d0<<"."<<endl;
    abort ();
  }
  else
    return table[i0];
}
