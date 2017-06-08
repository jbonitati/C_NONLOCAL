class vector_class<double> eig_vec (const double  eig_val,int dim,
				    const class array<double> &diag,
				    const class array<double> &off_diag,double precision)
{
  const int N = dim;

  const int Nm1 = N-1,Nm2 = N-2;

  class array<double>  sub_diag(N),super_diag(N),diag_shifted(N);
  class vector_class<double> X(N);

  class vector_class<double> vect_test(N);
  class vector_class<double> vect_test2(N);
  
  class vector_class<double> X_prev(N);
  

  X.random_vector ();
  X.normalization ();

  X_prev=X;
  
  double PRE=1e-5;

  int compt=0;

  do
    {  
      if (compt>100) { precision*=10;}

      for (unsigned int i = 1 ; i < N ; i++) sub_diag(i) = off_diag(i-1);
      for (unsigned int i = 0 ; i < N ; i++) super_diag(i) = off_diag(i),diag_shifted(i) = diag(i) - eig_val+PRE ;
      
      for (unsigned int i = 1 ; i < N ; i++)
	{
	  sub_diag(i) /= diag_shifted(i-1); 
	  diag_shifted(i) -= sub_diag(i)*super_diag(i-1);
	}
      
      for (unsigned int i = 1 ; i < N ; i++) X(i) -= sub_diag(i)*X(i-1);
      
      X(Nm1) /= diag_shifted(Nm1);
      for (unsigned int i = Nm2 ; i <= Nm2 ; i--) X(i) = (X(i) - super_diag(i)*X(i+1))/diag_shifted(i);
      
      X.normalization ();
      
      vect_test=X_prev-X;
      vect_test2=X_prev+X;
       
      
      X_prev=X;
    
      compt+=1;

    }
  
  while ( abs(vect_test.coef_max()) > precision  &&  abs(vect_test2.coef_max()) > precision   );
    
    
  return X;    
 
}







void Lanczos (int nb_vect_to_calc,double epsilon,int step_Lanczos
	      ,int dim_tot
	      ,class array <double> & eigen_val_final
	      ,class array <double> & vect_final
	      ,class matrix<double> & A)
{
 
  int dim_final_Lanczos;

  class vector_class<double> V_tmp (dim_tot); 
 
  class array<double> * eig_val; 
  eig_val=0;
  

  class vector_class<double> **eig_list; 
  eig_list=0;


  class array<double> * diag;
  class array<double> * off_diag;

  diag=0;
  off_diag=0;

  class array<double> previous_res(nb_vect_to_calc);
  class array<double> dif(nb_vect_to_calc);
   
  for (int i=0;i<=nb_vect_to_calc-1;i++)
    {
      previous_res(i)=0.0;
      dif(i)=0.0;
    }
	
  class vector_class<double> **V_Lanczos_tab; 
  V_Lanczos_tab=0;

   
  int i_previous;
  int i_lim;


  int dim;

  for (int i=0;i<=dim_tot-1;i++)
    {
      V_tmp(i)=0.0;
    }
  

  class vector_class<double> pivot(dim_tot);
  double norm_pivot=0.0;

  
  for (int i=0;i<=dim_tot-1;i++)
    {
      pivot(i)=drand48();
      norm_pivot+=pivot(i)*pivot(i);
    }

   
  for (int i=0;i<=dim_tot-1;i++)
    {
      pivot(i)*=1.0/sqrt(norm_pivot);
    }


  V_tmp=A*pivot;

 
  double  E1,E2;
      
  E1=pivot*V_tmp; 
      
  class vector_class <double>  vect2(dim_tot);
      
  vect2=V_tmp-E1*pivot;
  E2=vect2.dot_product_norm();
      
  vect2.normalization();
  
  diag=new class array<double>(dim_tot);
  off_diag=new class array<double>(dim_tot);
  
  
  (*diag)(0)=E1;
  (*off_diag)(0)=E2;
  
  
  dim=1;
   
	
  eig_val=new class array<double>(dim_tot);
	
  V_Lanczos_tab= new class vector_class<double> * [dim_tot];
  
  V_Lanczos_tab[0]=new class vector_class<double> (dim_tot);
  V_Lanczos_tab[1]=new class vector_class<double> (dim_tot);
  
  (*(V_Lanczos_tab[0]))=pivot;
  (*(V_Lanczos_tab[1]))=vect2;

  pivot=(*(V_Lanczos_tab[1]));
      
     


  i_previous=2;
  i_lim=(step_Lanczos<dim_tot)?  (step_Lanczos):(dim_tot)  ;

  
  bool keep_on;
  
  do 
    {
      keep_on=false;

      if ( i_previous>2 )
	{
	  for (int i=0;i<=dim-1;i++)
	    {
	      delete  eig_list[i];
	      eig_list[i]=0;
	    }
	  
	  delete [] eig_list;
	  eig_list=0;
	}
	 
	  
      for (int i=i_previous;i<=i_lim;i++)
	{	    

	   
	  for (int a=0;a<=dim_tot-1;a++)
	    {
	      V_tmp(a)=0.0;
	    }
	   	 
	      
	  V_tmp=A*pivot;


	  (*diag)(i-1)=(*(V_Lanczos_tab[i-1]))*V_tmp;
	      
	  if (i<=dim_tot-1)
	    {
	      V_Lanczos_tab[i]=new class vector_class<double> (dim_tot);
		  
	      (*(V_Lanczos_tab[i]))=V_tmp
		-(*diag)(i-1)*(*(V_Lanczos_tab[i-1]))-(*off_diag)(i-2)*(*(V_Lanczos_tab[i-2]));
	      
	      (*off_diag)(i-1)=(*(V_Lanczos_tab[i])).dot_product_norm();
		  
	      (*(V_Lanczos_tab[i])).normalization();
		  
	      // Graamm Schmidt
	      
	      for (int k=0;k<=i-1;k++)
		{
		  (*(V_Lanczos_tab[i]))-=((*(V_Lanczos_tab[i]))*(*(V_Lanczos_tab[k])))*(*(V_Lanczos_tab[k]));
		}
	      
	      (*(V_Lanczos_tab[i])).normalization();
		  
	      pivot=(*(V_Lanczos_tab[i]));
	    }		       
	      
	  dim+=1;
	  
	}
	  
	 
      class array<double> diag_small(dim),off_diag_small(dim),E_tab_small(dim);
	      
      for (unsigned int j = 0 ; j<=dim-1 ; j++) 
	{
	  diag_small(j) = (*diag)(j);
	  off_diag_small(j) = (*off_diag)(j);
	}
      
      off_diag_small(dim-1)=0.0;
	      
      //calcul des valeurs propres
	      
      Householder_QL2::all_eigenvalues_Lanc (dim,diag_small,off_diag_small,E_tab_small);
	  
      	  
      //iteration inverse pour les vecteurs
	  
      class vector_class<double> **eigenvectors;
	  
	  
      eigenvectors= new class vector_class<double> * [dim];
      
      for ( int z=0 ; z <=dim-1; z++)
	{
	  eigenvectors[z]= new class vector_class<double> (dim);
	      
	  *(eigenvectors[z])= eig_vec(E_tab_small(z),dim,diag_small,off_diag_small,1e-10); 
	}
	  
      	  
      eig_list= new class vector_class<double> * [dim];
	  
      for ( int z=0 ; z <=dim-1; z++)
	{ 
	  (*eig_val)(z)=E_tab_small(z);
	  eig_list[z] = new class vector_class<double> (dim_tot);
	      
	  class vector_class<double> &eig_vec = *(eig_list[z]);
	      
	  for ( int h=0 ; h<=dim_tot-1; h++)
	    {
	      eig_vec(h)=0.0;
	    }
	  
	  for ( int h=0 ; h<=dim-1; h++)
	    {
	      eig_vec+= (*(eigenvectors[z]))(h)*(*(V_Lanczos_tab[h]));
	    }
	}
      
      for (int j=0;j<=dim-1;j++)
	{
	  delete eigenvectors[j];
	  eigenvectors[j]=0;
	}
      
      delete [] eigenvectors;
      
      eigenvectors=0;
      
      cout<<"iteration "<<dim<<"/"<<dim_tot<<":"<<endl;
	   
	
      for (int k=0;k<=min(nb_vect_to_calc,dim)-1;k++)
	{
	  cout<<"eigenvalues Lanczos "<<E_tab_small(k)<<endl;
	  
	  eigen_val_final(k)=E_tab_small(k);
	        
	  dif(k)=E_tab_small(k)-previous_res(k);
	      
 
	  if ( abs(dif(k))>epsilon  )
	    {
	      keep_on=true;
	    }

	
      
	  previous_res(k)=E_tab_small(k);
	}
      
      i_previous=i_lim+1;
      i_lim=(i_previous+step_Lanczos<dim_tot)?  (i_previous+step_Lanczos):(dim_tot)  ;
      
      if (i_lim<i_previous) {i_previous=dim_tot+5;}  // pour arreter la boucle
            
      dim_final_Lanczos=dim;
	
    }

  
  while ( (i_previous<=dim_tot && keep_on==true) ||  (dim<nb_vect_to_calc) );


  for ( int i=0;i<=nb_vect_to_calc-1; i++)
    {
      double fact=1.0;
      
      if ( (*(eig_list[i]))(0) <1e-9  )
	{
	  fact=1.0;
	}

    for ( int h=0;h<=dim_tot-1; h++)
      {
        vect_final(i,h)=fact*(*(eig_list[i]))(h) ;
      }
    }



  for (int i=0;i<=dim-1;i++)
    {
      delete  eig_list[i];
      eig_list[i]=0;
    }
  
  delete [] eig_list;
  
  eig_list=0;
  
  
  for (int i=0;i<=dim_final_Lanczos-1;i++)
    {
      delete  V_Lanczos_tab[i];
      V_Lanczos_tab[i]=0;
    }
  
  delete [] V_Lanczos_tab;
  
  V_Lanczos_tab=0;
  
  
  delete diag;
  diag=0;
  
  delete off_diag;
  off_diag=0;
  
  
  delete eig_val;
  eig_val=0;
  
  	
}


class vector_class<complex < double> >  eig_vec_complex (const complex < double >  eig_val,int dim,
					    const class array<complex < double> >  &diag,
					    const class array<complex < double> >  &off_diag,double precision)
{
  const int N = dim;

  const int Nm1 = N-1,Nm2 = N-2;

  class array<complex < double> >   sub_diag(N),super_diag(N),diag_shifted(N);
  class vector_class<complex < double> >  X(N);

  class vector_class<complex < double> >  vect_test(N);
  class vector_class<complex < double> >  vect_test2(N);
  
  class vector_class<complex < double> >  X_prev(N);
  

  X.random_vector ();
  X.normalization ();

  X_prev=X;
  
  double PRE=1e-5;

  int compt=0;

  do
    {  
      if (compt>100) { precision*=10;}

      for (unsigned int i = 1 ; i < N ; i++) sub_diag(i) = off_diag(i-1);
      for (unsigned int i = 0 ; i < N ; i++) super_diag(i) = off_diag(i),diag_shifted(i) = diag(i) - eig_val+PRE ;
      
      for (unsigned int i = 1 ; i < N ; i++)
	{
	  sub_diag(i) /= diag_shifted(i-1); 
	  diag_shifted(i) -= sub_diag(i)*super_diag(i-1);
	}
      
      for (unsigned int i = 1 ; i < N ; i++) X(i) -= sub_diag(i)*X(i-1);
      
      X(Nm1) /= diag_shifted(Nm1);
      for (unsigned int i = Nm2 ; i <= Nm2 ; i--) X(i) = (X(i) - super_diag(i)*X(i+1))/diag_shifted(i);
      
      X.normalization ();
      
      vect_test=X_prev-X;
      vect_test2=X_prev+X;
       
      
      X_prev=X;
    
      compt+=1;

    }
  
  while ( abs(vect_test.coef_max()) > precision  &&  abs(vect_test2.coef_max()) > precision   );
    
    
  return X;    
 
}





void Lanczos_complex (int nb_vect_to_calc,double epsilon,int step_Lanczos
		      ,int dim_tot
		      ,class array <complex < double> >  & eigen_val_final
		      ,class array <complex < double> >  & vect_final
		      ,class matrix<complex < double> >  & A)
{
  
  int dim_final_Lanczos;

  class vector_class<complex < double> >  V_tmp (dim_tot); 
 
  class array<complex < double> >  * eig_val; 
  eig_val=0;
  

  class vector_class<complex < double> >  **eig_list; 
  eig_list=0;


  class array<complex < double> >  * diag;
  class array<complex < double> >  * off_diag;

  diag=0;
  off_diag=0;

  class array<complex < double> >  previous_res(nb_vect_to_calc);
  class array<complex < double> >  dif(nb_vect_to_calc);
   
  for (int i=0;i<=nb_vect_to_calc-1;i++)
    {
      previous_res(i)=0.0;
      dif(i)=0.0;
    }
	
  class vector_class<complex < double> >  **V_Lanczos_tab; 
  V_Lanczos_tab=0;

   
  int i_previous;
  int i_lim;


  int dim;

  for (int i=0;i<=dim_tot-1;i++)
    {
      V_tmp(i)=0.0;
    }
  

  class vector_class<complex < double> >  pivot(dim_tot);
 complex< double > norm_pivot=0.0;

  
  for (int i=0;i<=dim_tot-1;i++)
    {
      pivot(i)=drand48();
      norm_pivot+=pivot(i)*pivot(i);
    }

   
  for (int i=0;i<=dim_tot-1;i++)
    {
      pivot(i)*=1.0/sqrt(norm_pivot);
    }


  V_tmp=A*pivot;

 
 complex< double>  E1,E2;
      
  E1=pivot*V_tmp; 
      
  class vector_class <complex < double> >   vect2(dim_tot);
      
  vect2=V_tmp-E1*pivot;
  E2=vect2.dot_product_norm();
      
  vect2.normalization();
  
  diag=new class array<complex < double> > (dim_tot);
  off_diag=new class array<complex < double> > (dim_tot);
  
  
  (*diag)(0)=E1;
  (*off_diag)(0)=E2;
  
  
  dim=1;
   
	
  eig_val=new class array<complex < double> > (dim_tot);
	
  V_Lanczos_tab= new class vector_class<complex < double> >  * [dim_tot];
  
  V_Lanczos_tab[0]=new class vector_class<complex < double> >  (dim_tot);
  V_Lanczos_tab[1]=new class vector_class<complex < double> >  (dim_tot);
  
  (*(V_Lanczos_tab[0]))=pivot;
  (*(V_Lanczos_tab[1]))=vect2;

  pivot=(*(V_Lanczos_tab[1]));
      
     


  i_previous=2;
  i_lim=(step_Lanczos<dim_tot)?  (step_Lanczos):(dim_tot)  ;

  
  bool keep_on;
  
  do 
    {
      keep_on=false;

      if ( i_previous>2 )
	{
	  for (int i=0;i<=dim-1;i++)
	    {
	      delete  eig_list[i];
	      eig_list[i]=0;
	    }
	  
	  delete [] eig_list;
	  eig_list=0;
	}
	 
	  
      for (int i=i_previous;i<=i_lim;i++)
	{	    

	   
	  for (int a=0;a<=dim_tot-1;a++)
	    {
	      V_tmp(a)=0.0;
	    }
	   	 
	      
	  V_tmp=A*pivot;


	  (*diag)(i-1)=(*(V_Lanczos_tab[i-1]))*V_tmp;
	      
	  if (i<=dim_tot-1)
	    {
	      V_Lanczos_tab[i]=new class vector_class<complex < double> >  (dim_tot);
		  
	      (*(V_Lanczos_tab[i]))=V_tmp
		-(*diag)(i-1)*(*(V_Lanczos_tab[i-1]))-(*off_diag)(i-2)*(*(V_Lanczos_tab[i-2]));
	      
	      (*off_diag)(i-1)=(*(V_Lanczos_tab[i])).dot_product_norm();
		  
	      (*(V_Lanczos_tab[i])).normalization();
		  
	      // Graamm Schmidt
	      
	      for (int k=0;k<=i-1;k++)
		{
		  (*(V_Lanczos_tab[i]))-=((*(V_Lanczos_tab[i]))*(*(V_Lanczos_tab[k])))*(*(V_Lanczos_tab[k]));
		}
	      
	      (*(V_Lanczos_tab[i])).normalization();
		  
	      pivot=(*(V_Lanczos_tab[i]));
	    }		       
	      
	  dim+=1;
	  
	}
	  
	 
      class array<complex < double> >  diag_small(dim),off_diag_small(dim),E_tab_small(dim);
	      
      for (unsigned int j = 0 ; j<=dim-1 ; j++) 
	{
	  diag_small(j) = (*diag)(j);
	  off_diag_small(j) = (*off_diag)(j);
	}
      
      off_diag_small(dim-1)=0.0;
	      
      //calcul des valeurs propres
	      
      Householder_QL2::all_eigenvalues_Lanc (dim,diag_small,off_diag_small,E_tab_small);
	  
      	  
      //iteration inverse pour les vecteurs
	  
      class vector_class<complex < double> >  **eigenvectors;
	  
	  
      eigenvectors= new class vector_class<complex < double> >  * [dim];
      
      for ( int z=0 ; z <=dim-1; z++)
	{
	  eigenvectors[z]= new class vector_class<complex < double> >  (dim);
	      
	  *(eigenvectors[z])= eig_vec_complex(E_tab_small(z),dim,diag_small,off_diag_small,1e-10); 
	}
	  
      	  
      eig_list= new class vector_class<complex < double> >  * [dim];
	  
      for ( int z=0 ; z <=dim-1; z++)
	{ 
	  (*eig_val)(z)=E_tab_small(z);
	  eig_list[z] = new class vector_class<complex < double> >  (dim_tot);
	      
	  class vector_class<complex < double> >  &eig_vec = *(eig_list[z]);
	      
	  for ( int h=0 ; h<=dim_tot-1; h++)
	    {
	      eig_vec(h)=0.0;
	    }
	  
	  for ( int h=0 ; h<=dim-1; h++)
	    {
	      eig_vec+= (*(eigenvectors[z]))(h)*(*(V_Lanczos_tab[h]));
	    }
	}
      
      for (int j=0;j<=dim-1;j++)
	{
	  delete eigenvectors[j];
	  eigenvectors[j]=0;
	}
      
      delete [] eigenvectors;
      
      eigenvectors=0;
      
      cout<<"iteration "<<dim<<"/"<<dim_tot<<":"<<endl;
	   
	
      for (int k=0;k<=min(nb_vect_to_calc,dim)-1;k++)
	{
	  cout<<"eigenvalues Lanczos "<<E_tab_small(k)<<endl;
	  
	  eigen_val_final(k)=E_tab_small(k);
	        
	  dif(k)=E_tab_small(k)-previous_res(k);
	      
 
	  if ( abs(dif(k))>epsilon  )
	    {
	      keep_on=true;
	    }

	
      
	  previous_res(k)=E_tab_small(k);
	}
      
      i_previous=i_lim+1;
      i_lim=(i_previous+step_Lanczos<dim_tot)?  (i_previous+step_Lanczos):(dim_tot)  ;
      
      if (i_lim<i_previous) {i_previous=dim_tot+5;}  // pour arreter la boucle
            
      dim_final_Lanczos=dim;
	
    }

  
  while ( (i_previous<=dim_tot && keep_on==true) ||  (dim<nb_vect_to_calc) );


  for ( int i=0;i<=nb_vect_to_calc-1; i++)
    {
      double fact=1.0;
      
      if ( real( (*(eig_list[i]))(0) ) <1e-9  )
	{
	  fact=1.0;
	}

    for ( int h=0;h<=dim_tot-1; h++)
      {
        vect_final(i,h)=fact*(*(eig_list[i]))(h) ;
      }
    }



  for (int i=0;i<=dim-1;i++)
    {
      delete  eig_list[i];
      eig_list[i]=0;
    }
  
  delete [] eig_list;
  
  eig_list=0;
  
  
  for (int i=0;i<=dim_final_Lanczos-1;i++)
    {
      delete  V_Lanczos_tab[i];
      V_Lanczos_tab[i]=0;
    }
  
  delete [] V_Lanczos_tab;
  
  V_Lanczos_tab=0;
  
  
  delete diag;
  diag=0;
  
  delete off_diag;
  off_diag=0;
  
  
  delete eig_val;
  eig_val=0;
  
  	
}





void Lanczos_complex_with_pivot (int nb_vect_to_calc,double epsilon,int step_Lanczos
				 ,int dim_tot
				 ,class array <complex < double> >  & eigen_val_final
				 ,class array <complex < double> >  & vect_final
				 ,class matrix<complex < double> >  & A,class array <complex < double> >  & vect_pivot)
{
   
  class array <complex < double> >  overlap_tab(dim_tot);
   
  
 
  class array <complex < double> >   vect_pivot_original(dim_tot);

  for (int i=0;i<=dim_tot-1;i++)
    {
      vect_pivot_original(i)=vect_pivot(i);
    }


  int dim_final_Lanczos;

  class vector_class<complex < double> >  V_tmp (dim_tot); 
 
  class array<complex < double> >  * eig_val; 
  eig_val=0;
  

  class vector_class<complex < double> >  **eig_list; 
  eig_list=0;


  class array<complex < double> >  * diag;
  class array<complex < double> >  * off_diag;

  diag=0;
  off_diag=0;

  class array<complex < double> >  previous_res(nb_vect_to_calc);
  class array<complex < double> >  dif(nb_vect_to_calc);
   
  for (int i=0;i<=nb_vect_to_calc-1;i++)
    {
      previous_res(i)=0.0;
      dif(i)=0.0;
    }
	
  class vector_class<complex < double> >  **V_Lanczos_tab; 
  V_Lanczos_tab=0;

   
  int i_previous;
  int i_lim;


  int dim;

  for (int i=0;i<=dim_tot-1;i++)
    {
      V_tmp(i)=0.0;
    }
  

  class vector_class<complex < double> >  pivot(dim_tot);
 complex< double > norm_pivot=0.0;

  
  for (int i=0;i<=dim_tot-1;i++)
    {
      pivot(i)=vect_pivot(i);
      norm_pivot+=pivot(i)*pivot(i);
    }

   
  for (int i=0;i<=dim_tot-1;i++)
    {
      pivot(i)*=1.0/sqrt(norm_pivot);
    }


  V_tmp=A*pivot;

 
 complex< double>  E1,E2;
      
  E1=pivot*V_tmp; 
      
  class vector_class <complex < double> >   vect2(dim_tot);
      
  vect2=V_tmp-E1*pivot;
  E2=vect2.dot_product_norm();
      
  vect2.normalization();
  
  diag=new class array<complex < double> > (dim_tot);
  off_diag=new class array<complex < double> > (dim_tot);
  
  
  (*diag)(0)=E1;
  (*off_diag)(0)=E2;
  
  
  dim=1;
   
	
  eig_val=new class array<complex < double> > (dim_tot);
	
  V_Lanczos_tab= new class vector_class<complex < double> >  * [dim_tot];
  
  V_Lanczos_tab[0]=new class vector_class<complex < double> >  (dim_tot);
  V_Lanczos_tab[1]=new class vector_class<complex < double> >  (dim_tot);
  
  (*(V_Lanczos_tab[0]))=pivot;
  (*(V_Lanczos_tab[1]))=vect2;

  pivot=(*(V_Lanczos_tab[1]));
      
     


  i_previous=2;
  i_lim=(step_Lanczos<dim_tot)?  (step_Lanczos):(dim_tot)  ;

  
  bool keep_on;
  
  do 
    {
      keep_on=false;

      if ( i_previous>2 )
	{
	  for (int i=0;i<=dim-1;i++)
	    {
	      delete  eig_list[i];
	      eig_list[i]=0;
	    }
	  
	  delete [] eig_list;
	  eig_list=0;
	}
	 
	  
      for (int i=i_previous;i<=i_lim;i++)
	{	    

	   
	  for (int a=0;a<=dim_tot-1;a++)
	    {
	      V_tmp(a)=0.0;
	    }
	   	 
	      
	  V_tmp=A*pivot;


	  (*diag)(i-1)=(*(V_Lanczos_tab[i-1]))*V_tmp;
	      
	  if (i<=dim_tot-1)
	    {
	      V_Lanczos_tab[i]=new class vector_class<complex < double> >  (dim_tot);
		  
	      (*(V_Lanczos_tab[i]))=V_tmp
		-(*diag)(i-1)*(*(V_Lanczos_tab[i-1]))-(*off_diag)(i-2)*(*(V_Lanczos_tab[i-2]));
	      
	      (*off_diag)(i-1)=(*(V_Lanczos_tab[i])).dot_product_norm();
		  
	      (*(V_Lanczos_tab[i])).normalization();
		  
	      // Graamm Schmidt
	      
	      for (int k=0;k<=i-1;k++)
		{
		  (*(V_Lanczos_tab[i]))-=((*(V_Lanczos_tab[i]))*(*(V_Lanczos_tab[k])))*(*(V_Lanczos_tab[k]));
		}
	      
	      (*(V_Lanczos_tab[i])).normalization();
		  
	      pivot=(*(V_Lanczos_tab[i]));
	    }		       
	      
	  dim+=1;
	  
	}
	  
	 
      class array<complex < double> >  diag_small(dim),off_diag_small(dim),E_tab_small(dim);
	      
      for (unsigned int j = 0 ; j<=dim-1 ; j++) 
	{
	  diag_small(j) = (*diag)(j);
	  off_diag_small(j) = (*off_diag)(j);
	}
      
      off_diag_small(dim-1)=0.0;
	      
      //calcul des valeurs propres
	      
      Householder_QL2::all_eigenvalues_Lanc (dim,diag_small,off_diag_small,E_tab_small);
	  
      	  
      //iteration inverse pour les vecteurs
	  
      class vector_class<complex < double> >  **eigenvectors;
	  
	  
      eigenvectors= new class vector_class<complex < double> >  * [dim];
      
      for ( int z=0 ; z <=dim-1; z++)
	{
	  eigenvectors[z]= new class vector_class<complex < double> >  (dim);
	      
	  *(eigenvectors[z])= eig_vec_complex(E_tab_small(z),dim,diag_small,off_diag_small,1e-10); 
	}
	  
      	  
      eig_list= new class vector_class<complex < double> >  * [dim];
	  
      for ( int z=0 ; z <=dim-1; z++)
	{ 
	  (*eig_val)(z)=E_tab_small(z);
	  eig_list[z] = new class vector_class<complex < double> >  (dim_tot);
	      
	  class vector_class<complex < double> >  &eig_vec = *(eig_list[z]);
	      
	  for ( int h=0 ; h<=dim_tot-1; h++)
	    {
	      eig_vec(h)=0.0;
	    }
	  
	  for ( int h=0 ; h<=dim-1; h++)
	    {
	      eig_vec+= (*(eigenvectors[z]))(h)*(*(V_Lanczos_tab[h]));
	    }

	  
	  overlap_tab(z)=0.0;

	  for (int pl=0;pl<=dim_tot-1;pl++)
	    {
	      overlap_tab(z)+=eig_vec(pl)*vect_pivot_original(pl);
	    }


	  overlap_tab(z)=abs(overlap_tab(z));

	}
      
      for (int j=0;j<=dim-1;j++)
	{
	  delete eigenvectors[j];
	  eigenvectors[j]=0;
	}
      
      delete [] eigenvectors;
      
      eigenvectors=0;
      
      cout<<"iteration "<<dim<<"/"<<dim_tot<<":"<<endl;
	   
	
      for (int k=0;k<=min(nb_vect_to_calc,dim)-1;k++)
	{
	  cout<<"eigenvalues Lanczos "<<E_tab_small(k)<<" overlap "<<overlap_tab(k)<<endl;
	  
	  eigen_val_final(k)=E_tab_small(k);
	        
	  dif(k)=E_tab_small(k)-previous_res(k);
	      
 
	  if ( abs(dif(0))>epsilon  )
	    {
	      keep_on=true;
	    }

	
      
	  previous_res(k)=E_tab_small(k);
	}
      
      i_previous=i_lim+1;
      i_lim=(i_previous+step_Lanczos<dim_tot)?  (i_previous+step_Lanczos):(dim_tot)  ;
      
      if (i_lim<i_previous) {i_previous=dim_tot+5;}  // pour arreter la boucle
            
      dim_final_Lanczos=dim;
	
    }

  
  while ( (i_previous<=dim_tot && keep_on==true) ||  (dim<nb_vect_to_calc) );


  for ( int i=0;i<=nb_vect_to_calc-1; i++)
    {
      double fact=1.0;
      
      if ( real( (*(eig_list[i]))(0) ) <1e-9  )
	{
	  fact=1.0;
	}

    for ( int h=0;h<=dim_tot-1; h++)
      {
        vect_final(i,h)=fact*(*(eig_list[i]))(h) ;
      }
    }



  for (int i=0;i<=dim-1;i++)
    {
      delete  eig_list[i];
      eig_list[i]=0;
    }
  
  delete [] eig_list;
  
  eig_list=0;
  
  
  for (int i=0;i<=dim_final_Lanczos-1;i++)
    {
      delete  V_Lanczos_tab[i];
      V_Lanczos_tab[i]=0;
    }
  
  delete [] V_Lanczos_tab;
  
  V_Lanczos_tab=0;
  
  
  delete diag;
  diag=0;
  
  delete off_diag;
  off_diag=0;
  
  
  delete eig_val;
  eig_val=0;
  
  	
}




