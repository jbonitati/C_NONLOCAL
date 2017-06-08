// $Id: arrays.H,v 1.6 2002/12/14 17:06:04 lmichel Exp $

// templates for arrays of dimension 1,2,3,4,5,6,7.
// test of indices included.

struct pair_str
{
  int s0,s1;
};

template <class C>
class array
{
public:

  array (const unsigned int d0_c) 
    : d0(d0_c), d1(0), d2(0), d3(0), d4(0), d5(0), d6(0) ,d7(0),d8(0)
  {
    if (d0 == 0)
    {
      cout<<" d0 Dimension is zero."<<endl;
      abort();
    }
    

    table = new C [d0];
  }

  array (const unsigned int d0_c,const unsigned int d1_c) 
    : d0(d0_c), d1(d1_c), d2(0), d3(0), d4(0), d5(0), d6(0),d7(0) ,d8(0)
  {
    if ((d0 == 0) || (d1 == 0))
    {
      cout<<"Dimension is zero."<<" "<<d0<<" "<<d1<<endl;
      abort();
    }
    
    table = new C [d0*d1];
  }

  array (const unsigned int d0_c,const unsigned int d1_c,const unsigned int d2_c) 
    : d0(d0_c), d1(d1_c), d2(d2_c), d3(0), d4(0), d5(0), d6(0),d7(0) ,d8(0)
  { 
    if ((d0 == 0) || (d1 == 0) || (d2 == 0))
    {
      cout<<"Dimension is zero."<<" "<<d0<<" "<<d1<<" "<<d2<<endl;
      abort();
    }

    table = new C [d0*d1*d2];
  }
  
  array (const unsigned int d0_c,const unsigned int d1_c,const unsigned int d2_c,const unsigned int d3_c) 
    : d0(d0_c), d1(d1_c), d2(d2_c), d3(d3_c), d4(0), d5(0), d6(0),d7(0) ,d8(0)
  {
    if ((d0 == 0) || (d1 == 0) || (d2 == 0) || (d3 == 0))

    {
      cout<<"Dimension is zero."<<endl;
      abort();
    }

    table = new C [d0*d1*d2*d3];
  }

  array (const unsigned int d0_c,const unsigned int d1_c,const unsigned int d2_c,
         const unsigned int d3_c,const unsigned int d4_c) 
    : d0(d0_c), d1(d1_c), d2(d2_c), d3(d3_c), d4(d4_c), d5(0), d6(0),d7(0) ,d8(0)
  {
    if ((d0 == 0) || (d1 == 0) || (d2 == 0) || (d3 == 0) || (d4 == 0))
    {
      cout<<"Dimension is zero."<<endl;
      abort();
    }

    table = new C [d0*d1*d2*d3*d4];
  }
  
  array (const unsigned int d0_c,const unsigned int d1_c,const unsigned int d2_c,
         const unsigned int d3_c,const unsigned int d4_c,const unsigned int d5_c) 
    : d0(d0_c), d1(d1_c), d2(d2_c), d3(d3_c), d4(d4_c), d5(d5_c), d6(0),d7(0) ,d8(0)
  {
    if ((d0 == 0) || (d1 == 0) || (d2 == 0) || (d3 == 0) || (d4 == 0) || (d5 == 0))
    {
      cout<<"Dimension is zero."<<endl;
      abort();
    }

    table = new C [d0*d1*d2*d3*d4*d5];
  }

  array (const unsigned int d0_c,const unsigned int d1_c,const unsigned int d2_c,
         const unsigned int d3_c,const unsigned int d4_c,const unsigned int d5_c,const unsigned int d6_c) 
    : d0(d0_c), d1(d1_c), d2(d2_c), d3(d3_c), d4(d4_c), d5(d5_c), d6(d6_c),d7(0) ,d8(0)
  {
    if ((d0 == 0) || (d1 == 0) || (d2 == 0) || (d3 == 0) || (d4 == 0) || (d5 == 0) || (d6 == 0))
    {
      cout<<"Dimension is zero."<<endl;
      abort();
    }


    table = new C [d0*d1*d2*d3*d4*d5*d6];
  }

  array (const unsigned int d0_c,const unsigned int d1_c,const unsigned int d2_c,
         const unsigned int d3_c,const unsigned int d4_c,const unsigned int d5_c
	 ,const unsigned int d6_c,const unsigned int d7_c,const unsigned int d8_c)
 
    : d0(d0_c), d1(d1_c), d2(d2_c), d3(d3_c), d4(d4_c), d5(d5_c), d6(d6_c) ,d7(d7_c), d8(d8_c)
  {
    if ((d0 == 0) || (d1 == 0) || (d2 == 0) || (d3 == 0) || (d4 == 0) || (d5 == 0) || (d6 == 0)
	|| (d7 == 0) || (d8 == 0)  )
    {
      cout<<"Dimension is zero."<<endl;
      abort();
    }


    table = new C [d0*d1*d2*d3*d4*d5*d6*d7*d8];
  }


   array (const unsigned int d0_c,const unsigned int d1_c,const unsigned int d2_c,
         const unsigned int d3_c,const unsigned int d4_c,const unsigned int d5_c
         ,const unsigned int d6_c,const unsigned int d7_c)

    : d0(d0_c), d1(d1_c), d2(d2_c), d3(d3_c), d4(d4_c), d5(d5_c), d6(d6_c) ,d7(d7_c), d8(0)
  {
    if ((d0 == 0) || (d1 == 0) || (d2 == 0) || (d3 == 0) || (d4 == 0) || (d5 == 0) || (d6 == 0)
        || (d7 == 0) )
    {
      cout<<"Dimension is zero."<<endl;
      abort();
    }


    table = new C [d0*d1*d2*d3*d4*d5*d6*d7];
  }





  ~array (void)
  {
    delete [] table;
  }
  
  C *table;

  C & operator [] (const struct pair_str &pair) const;
  C & operator () (const int X,const struct pair_str &pair) const;


  C & operator [] (const unsigned int i0) const;

  C & operator () (const unsigned int i0) const;
  
  C & operator () (const unsigned int i0,const unsigned int i1) const;
  
  C & operator () (const unsigned int i0,const unsigned int i1,const unsigned int i2) const;

  C & operator () (const unsigned int i0,const unsigned int i1,
                   const unsigned int i2,const unsigned int i3) const;

  C & operator () (const unsigned int i0,const unsigned int i1,
                   const unsigned int i2,const unsigned int i3,const unsigned int i4) const;

  C & operator () (const unsigned int i0,const unsigned int i1,
                   const unsigned int i2,const unsigned int i3,
                   const unsigned int i4,const unsigned int i5) const;

  C & operator () (const unsigned int i0,const unsigned int i1,
                   const unsigned int i2,const unsigned int i3,
                   const unsigned int i4,const unsigned int i5,const unsigned int i6) const;

 C & operator () (const unsigned int i0,const unsigned int i1,
                   const unsigned int i2,const unsigned int i3,
                  const unsigned int i4,const unsigned int i5,const unsigned int i6
                  ,const unsigned int i7) const;

 
 C & operator () (const unsigned int i0,const unsigned int i1,
                   const unsigned int i2,const unsigned int i3,
		  const unsigned int i4,const unsigned int i5,const unsigned int i6
		  ,const unsigned int i7,const unsigned int i8) const;


  const unsigned int d0,d1,d2,d3,d4,d5,d6,d7,d8;
};

template <class C>
C & array<C>::operator [] (const struct pair_str &pair) const

{
  if (d2 > 0)
  {
    cout<<"1 component missing."<<endl;
    abort ();
  }
  else if (pair.s0 >= d0)
  {
    cout<<"Problem with pair.s0="<<pair.s0<<" >= d0="<<d0<<"."<<endl;
    abort ();
  }
  else if (pair.s1 >= d1)
  {
    cout<<"Problem with pair.s1="<<pair.s1<<" >= d1="<<d1<<"."<<endl;
    abort ();
  }
  else
    return table[pair.s1 + d1*pair.s0];
} 

template <class C>

C & array<C>::operator () (const int X,const struct pair_str &pair) const
{
  if (d2 == 0)
  {
    cout<<"1 component too much."<<endl;
    abort ();
  }
  else if (X >= d0)
  {
    cout<<"Problem with X="<<X<<" >= d0="<<d0<<"."<<endl;
    abort ();
  }
  else if (pair.s0 >= d1)
  {
    cout<<"Problem with pair.s0="<<pair.s0<<" >= d1="<<d1<<"."<<endl;
    abort ();
  }
  else if (pair.s1 >= d2)
  {
    cout<<"Problem with pair.s1="<<pair.s1<<" >= d2="<<d2<<"."<<endl;
    abort ();
  }
  else
    return table[pair.s1 + d2*(pair.s0 + d1*X)];
} 




template <class C>
C & array<C>::operator [] (const unsigned int i0) const
{
  if (d6 > 0)
  {
    cout<<"6 components missing. d6= "<<d6<<endl;

    cout<<"dimension= d0= "<<d0<<endl;

    abort ();
  }
  else if (d5 > 0)
  {
    cout<<"5 components missing."<<endl;
    abort ();
  }
  else if (d4 > 0)
  {
    cout<<"4 components missing."<<endl;
    abort ();
  }
  else if (d3 > 0)
  {
    cout<<"3 components missing."<<endl;
    abort ();
  } 
  else if (d2 > 0)
  {
    cout<<"2 components missing."<<endl;
    abort ();
  }
  else if (d1 > 0)
  {
    cout<<"1 component missing."<<endl;
    abort ();
  }
  else if (i0 >= d0)
  {
    cout<<"Problem with i0="<<i0<<" >= d0="<<d0<<"."<<endl;
   
    abort ();
  }
  else
    return table[i0];
}

template <class C>
C & array<C>::operator () (const unsigned int i0) const
{
  if (d6 > 0)
  {
    
    cout<<"6 components missing. d6= "<<d6<<endl;

    cout<<"dimension= d0= "<<d0<<endl;

    abort ();
  }

  else if (d5 > 0)
  {
    cout<<"5 components missing."<<endl;
    abort ();
  }
  else if (d4 > 0)
  {
    cout<<"4 components missing."<<endl;
    abort ();
  }
  else if (d3 > 0)
  {
    cout<<"3 components missing."<<endl;
    abort ();
  } 
  else if (d2 > 0)
  {
    cout<<"2 components missing."<<endl;
    abort ();
  }
  else if (d1 > 0)
  {
    cout<<"1 component missing."<<endl;
    abort ();
  }
  else if (i0 >= d0)
  {
    cout<<"Problem with i0="<<i0<<" >= d0="<<d0<<"."<<endl;
    
    cout<<"error b "<<endl;
    

    abort ();
  }
  else
    return table[i0];
}

template <class C>
C & array<C>::operator () (const unsigned int i0,const unsigned int i1) const
{
  if (d6 > 0)
  {
    cout<<"5 components missing."<<endl;
    abort ();
  }
  else if (d5 > 0)
  {
    cout<<"4 components missing."<<endl;
    abort ();
  }
  else if (d4 > 0)
  {
    cout<<"3 components missing."<<endl;
    abort ();
  }
  else if (d3 > 0)
  {
    cout<<"2 components missing."<<endl;
    abort ();
  } 
  else if (d2 > 0)
  {
    cout<<"1 component missing."<<endl;
    abort ();
  } 
  else if (d1 == 0)
    {
      cout<<"1 component too much."<<endl;
      abort ();
    }
  else if (i0 >= d0)
  {
    cout<<"Problem with i0="<<i0<<" >= d0="<<d0<<"."<<endl;
    
    cout<<"error c "<<endl;

    cout<<"error dims are "<<d0<<" "<<d1<<endl;

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

template <class C>
C & array<C>::operator () (const unsigned int i0,const unsigned int i1,const unsigned int i2) const
{
  if (d6 > 0)
  {
    cout<<"4 components missing."<<endl;
    abort ();
  }
  else if (d5 > 0)
  {
    cout<<"3 components missing."<<endl;
    abort ();
  }
  else if (d4 > 0)
  {
    cout<<"2 components missing."<<endl;
    abort ();
  }
  else if (d3 > 0)
  {
    cout<<"1 component missing."<<endl;
    abort ();
  }
  else if (d1 == 0)
  {
    cout<<"2 components too many."<<endl;
    abort ();
  } 
  else if (d2 == 0)
  {
    cout<<"1 component too much."<<endl;
    abort ();
  } 
  else if (i0 >= d0)
  {
    cout<<"Problem with i0="<<i0<<" >= d0="<<d0<<"."<<endl;
    
    cout<<"error d "<<endl;

    abort ();
  }
  else if (i1 >= d1)
  {
    cout<<"Problem with i1="<<i1<<" >= d1="<<d1<<"."<<endl;
    abort ();
  }
  else if (i2 >= d2)
  {
    cout<<"Problem with i2="<<i2<<" >= d2="<<d2<<"."<<endl;
    abort ();
  }
  else
    return table[i2 + d2*(i1 + d1*i0)];
}

template <class C>
C & array<C>::operator () (const unsigned int i0,const unsigned int i1,
                           const unsigned int i2,const unsigned int i3) const
{
  if (d6 > 0)
  {
    cout<<"3 components missing."<<endl;
    abort ();
  }
  else if (d5 > 0)
  {
    cout<<"2 components missing."<<endl;
    abort ();
  }
  else if (d4 > 0)
  {
    cout<<"1 component missing."<<endl;
    abort ();
  } 
  else if (d1 == 0)
  {
    cout<<"3 components too many."<<endl;
    abort ();
  } 
  else if (d2 == 0)
  {
    cout<<"2 components too many."<<endl;
    abort ();
  } 
  else if (d3 == 0)
  {
    cout<<"1 component too much."<<endl;
    abort ();
  }
  else if (i0 >= d0)
  {
    cout<<"Problem with i0="<<i0<<" >= d0="<<d0<<"."<<endl;
   
    cout<<"error e "<<endl;

    abort ();
  }
  else if (i1 >= d1)
  {
    cout<<"Problem with i1="<<i1<<" >= d1="<<d1<<"."<<endl;
    abort ();
  }
  else if (i2 >= d2)
  {
    cout<<"Problem with i2="<<i2<<" >= d2="<<d2<<"."<<endl;

    abort ();
  }
  else if (i3 >= d3)
  {
    cout<<"Problem with i3="<<i3<<" >= d3="<<d3<<"."<<endl;
    abort ();
  }
  else
    return table[i3 + d3*(i2 + d2*(i1 + d1*i0))];
}

template <class C>
C & array<C>::operator () (const unsigned int i0,const unsigned int i1,
                           const unsigned int i2,const unsigned int i3,const unsigned int i4) const
{
  if (d6 > 0)
  {
    cout<<"2 components missing."<<endl;
    abort ();
  }

  else if (d5 > 0)
  {
    cout<<"1 components missing."<<endl;
    abort ();
  }
  else if (d1 == 0)
  {
    cout<<"4 components too many."<<endl;
    abort ();
  }
  else if (d2 == 0)
  {
    cout<<"3 components too many."<<endl;
    abort ();
  } 
  else if (d3 == 0)
  {
    cout<<"2 components too many."<<endl;
    abort ();
  } 
  else if (d4 == 0)
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
  else if (i2 >= d2)
  {
    cout<<"Problem with i2="<<i2<<" >= d2="<<d2<<"."<<endl;
    abort ();
  }
  else if (i3 >= d3)
  {
    cout<<"Problem with i3="<<i3<<" >= d3="<<d3<<"."<<endl;
    abort ();
  }
  else if (i4 >= d4)
  {
    cout<<"Problem with i4="<<i4<<" >= d4="<<d4<<"."<<endl;
    abort ();
  }
  else
    return table[i4 + d4*(i3 + d3*(i2 + d2*(i1 + d1*i0)))];
}

template <class C>
C & array<C>::operator () (const unsigned int i0,const unsigned int i1,
                           const unsigned int i2,const unsigned int i3,
                           const unsigned int i4,const unsigned int i5) const
{
  if (d6 > 0)
  {
    cout<<"1 component missing."<<endl;
    abort ();
  }
  else if (d1 == 0)
  {
    cout<<"5 components too many."<<endl;
    abort ();
  }
  else if (d2 == 0)
  {
    cout<<"4 components too many."<<endl;
    abort ();
  } 
  else if (d3 == 0)
  {
    cout<<"3 components too many."<<endl;
    abort ();
  } 
  else if (d4 == 0)
  {
    cout<<"2 components too many."<<endl;
    abort ();
  }
  else if (d5 == 0)
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
  else if (i2 >= d2)
  {
    cout<<"Problem with i2="<<i2<<" >= d2="<<d2<<"."<<endl;
    abort ();
  }

  else if (i3 >= d3)
  {
    cout<<"Problem with i3="<<i3<<" >= d3="<<d3<<"."<<endl;
    abort ();
  }
  else if (i4 >= d4)
  {

    cout<<"Problem with i4="<<i4<<" >= d4="<<d4<<"."<<endl;
    abort ();
  }
  else if (i5 >= d5)
  {
    cout<<"Problem with i5="<<i5<<" >= d5="<<d5<<"."<<endl;
    abort ();
  }
  else
    return table[i5 + d5*(i4 + d4*(i3 + d3*(i2 + d2*(i1 + d1*i0))))];
}


template <class C>
C & array<C>::operator () (const unsigned int i0,const unsigned int i1,
                           const unsigned int i2,const unsigned int i3,
                           const unsigned int i4,const unsigned int i5,const unsigned int i6) const
{
  if (d1 == 0)
  {
    cout<<"6 components too many."<<endl;
    abort ();
  }
  else if (d2 == 0)
  {
    cout<<"5 components too many."<<endl;
    abort ();
  } 
  else if (d3 == 0)
  {
    cout<<"4 components too many."<<endl;
    abort ();
  } 
  else if (d4 == 0)
  {
    cout<<"3 components too many."<<endl;
    abort ();
  }
  else if (d5 == 0)
  {
    cout<<"2 components too many."<<endl;
    abort ();
  }

  else if (d6 == 0)
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
  else if (i2 >= d2)
  {
    cout<<"Problem with i2="<<i2<<" >= d2="<<d2<<"."<<endl;
    abort ();
  }
  else if (i3 >= d3)
  {
    cout<<"Problem with i3="<<i3<<" >= d3="<<d3<<"."<<endl;
    abort ();
  }
  else if (i4 >= d4)
  {
    cout<<"Problem with i4="<<i4<<" >= d4="<<d4<<"."<<endl;
    abort ();
  }
  else if (i5 >= d5)
  {
    cout<<"Problem with i5="<<i5<<" >= d5="<<d5<<"."<<endl;
    abort ();
  }
  else if (i6 >= d6)
  {
    cout<<"Problem with i6="<<i6<<" >= d6="<<d6<<"."<<endl;
    abort ();
  }
  else
    return table[i6 + d6*(i5 + d5*(i4 + d4*(i3 + d3*(i2 + d2*(i1 + d1*i0)))))];
}

template <class C>
C & array<C>::operator () (const unsigned int i0,const unsigned int i1,
                           const unsigned int i2,const unsigned int i3,
                           const unsigned int i4,const unsigned int i5,const unsigned int i6
                           ,const unsigned int i7) const
{
  if (d1 == 0)
  {
    cout<<"6 components too many."<<endl;
    abort ();
  }

   else if (d2 == 0)
  {
    cout<<"5 components too many."<<endl;
    abort ();
  }
  else if (d3 == 0)
  {
    cout<<"4 components too many."<<endl;
    abort ();
  }
  else if (d4 == 0)
  {
    cout<<"3 components too many."<<endl;
    abort ();
  }
  
   else if (d5 == 0)
  {
    cout<<"2 components too many."<<endl;
    abort ();
  }

  else if (d6 == 0)
  {
    cout<<"1 component too much."<<endl;
    abort ();
  }

 else if (d7 == 0)
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
  else if (i2 >= d2)
  {
    cout<<"Problem with i2="<<i2<<" >= d2="<<d2<<"."<<endl;
    abort ();
  }
  else if (i3 >= d3)
  {
    cout<<"Problem with i3="<<i3<<" >= d3="<<d3<<"."<<endl;
    abort ();
  }
  else if (i4 >= d4)
  {
    cout<<"Problem with i4="<<i4<<" >= d4="<<d4<<"."<<endl;
    abort ();
  }

 else if (i5 >= d5)
  {
    cout<<"Problem with i5="<<i5<<" >= d5="<<d5<<"."<<endl;
    abort ();
  }

  else if (i6 >= d6)
  {
    cout<<"Problem with i6="<<i6<<" >= d6="<<d6<<"."<<endl;
    abort ();
  }

 else if (i7 >= d7)
  {
    cout<<"Problem with i7="<<i7<<" >= d7="<<d7<<"."<<endl;
    abort ();
  }

 else
    return table[i7+ d7*( i6 + d6*(i5 + d5*(i4 + d4*(i3 + d3*(i2 + d2*(i1 + d1*i0))))))];


}

template <class C>
C & array<C>::operator () (const unsigned int i0,const unsigned int i1,
                           const unsigned int i2,const unsigned int i3,
                           const unsigned int i4,const unsigned int i5,const unsigned int i6
			   ,const unsigned int i7,const unsigned int i8) const
{
  if (d1 == 0)
  {
    cout<<"6 components too many."<<endl;
    abort ();
  }
  else if (d2 == 0)
  {
    cout<<"5 components too many."<<endl;
    abort ();
  } 
  else if (d3 == 0)
  {
    cout<<"4 components too many."<<endl;
    abort ();
  } 
  else if (d4 == 0)
  {
    cout<<"3 components too many."<<endl;
    abort ();
  }
  else if (d5 == 0)
  {
    cout<<"2 components too many."<<endl;
    abort ();
  }

  else if (d6 == 0)
  {
    cout<<"1 component too much."<<endl;
    abort ();
  }

 else if (d7 == 0)
  {
    cout<<"1 component too much."<<endl;
    abort ();
  }

 else if (d8 == 0)
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
  else if (i2 >= d2)
  {
    cout<<"Problem with i2="<<i2<<" >= d2="<<d2<<"."<<endl;
    abort ();
  }
  else if (i3 >= d3)
  {
    cout<<"Problem with i3="<<i3<<" >= d3="<<d3<<"."<<endl;
    abort ();
  }
  else if (i4 >= d4)
  {
    cout<<"Problem with i4="<<i4<<" >= d4="<<d4<<"."<<endl;
    abort ();
  }
  else if (i5 >= d5)
  {
    cout<<"Problem with i5="<<i5<<" >= d5="<<d5<<"."<<endl;
    abort ();
  }
  
  else if (i6 >= d6)
  {
    cout<<"Problem with i6="<<i6<<" >= d6="<<d6<<"."<<endl;
    abort ();
  }

 else if (i7 >= d7)
  {
    cout<<"Problem with i7="<<i7<<" >= d7="<<d7<<"."<<endl;
    abort ();
  }


  else if (i8 >= d8)
    {
    cout<<"Problem with i8="<<i8<<" >= d8="<<d8<<"."<<endl;
    abort ();
  }
  

  else
    return table[i8+d8*(i7+ d7*( i6 + d6*(i5 + d5*(i4 + d4*(i3 + d3*(i2 + d2*(i1 + d1*i0)))))))];
}

