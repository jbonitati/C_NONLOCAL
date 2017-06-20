#include <cmath>
#include <gsl/gsl_sf.h>

const double hbarc = 197.3269718;
const double mass_unit=931.494;
const double E2HC = 0.00729927;
const arma::cx_double I(0.0,1.0);

/*Returns the Coulomb potential for particles at distance R*/
double coulomb_potential(double R)
{
  double rc = 1.25; //don't know how to calculate this
  double Rcoul = rc*pow((targ.getM()),1.0/3.0);
  double z1 = targ.getZ();
  double z2 = proj.getZ();

  if (R<Rcoul)
  {
    return ((z1*z2*1.43997)/(2.0))*(3.0-pow(R,2.0)/(pow(Rcoul,2)));
  }
  else
  {
    return (z1*z2*1.43997)/R;
  }
}  

//calculates the coupling potential between two channels
//****incomplete*****
double couplingPotential(double r1, double r2, Channel * c1, Channel * c2){
  if(c1 == c2){
    return 0;
  }else{
    return beta * (pot.totalPotential(r1, targ, c1) 
      + nlpot.totalPotential(r1,r2,targ,c1));
  }
}

arma::mat<arma::cx_mat> System::cmatrixCalc(){
  //the C matrix is made up of a different matrix for each pair of channels
  arma::mat<arma::cx_mat> cmatrix(num_chanels,num_channels);
  int N = basis_size;
  
  //to iterate through the channels vector:
  std::vector<Channel *>::iterator c1, c2;
  int c = 0, cprime = 0;
  
  //TLmatrix uses the T_0, L(0), and energy operators
  //Vmatrix uses the potential, centrifugal, and channel energy operators
  for(c1 = channels.begin(); c1 != channels.end(); c1++){
    
    arma::cx_mat TLmatrix(N, N);
    double mu = (*c1)->getMu();
    double TLcoeff = hbarc*hbarc/(2*mu);
    
    for(int i=0;i<N;i++){ 
      double xi = lagrangeBasis.x(i);
      double ri = a*xi;
      
      TLmatrix(i,i) = TLcoeff*((4.*N*N+4.*N+3.0)*xi*(1.-xi)-6.*xi+1.)/
        (3.*ri*ri*((1.-xi)*(1.-xi))) - energy;
        
      for(int j=i+1;j<N;j++){
        double xj = lagrangeBasis.x(j);
        double rj = a*xj;
        
        double part3=pow(-1.,i+j)/(a*a*sqrt(xi*xj*(1.-xi)*(1.-xj)));
        double part4=(N*N*1.+N*1.0+1.0+(xi+xj-2.*xi*xj)/
          ((xi-xj)*(xi-xj))-1./(1.-xi)-1./(1.-xj));
          
        TLmatrix(i,j)=TLcoeff*part3*part4;
        TLmatrix(j,i)=TLmatrix(i,j);
      }
    } 
    
    for(c2 = channels.begin(); c2 != channels.end; c2++){
      arma::cx_mat Vmatrix(N,N, arma::fill::zeros);
      if((*c1)==(*c2)){
        //include TLmatrix as well as Vmatrix
        
        for(int i = 0; i < N; i++){
          double xi = lagrangeBasis.x(i);
          double ri = a*xi;
          
          //the diagonal terms include the central potential, energy 
          //eigenvalue, and local potential
          Vmatrix(i,i) = (*c1)->central_potential(ri,proj,targ) 
            + (*c1)->getE()
            + pot.totalPotential(ri, targ, (*c1))
            + coulomb_potential(ri);
        
          for(int j = 0; j < N; j++){
            double xj = lagrangeBasis.x(j);
            double rj = a*xj;
            
            Vmatrix(i,j) += a*sqrt(lagrangeBasis.w(i)*lagrangeBasis.w(j)) 
              * nlpot.totalPotential(ri, ri, targ, (*c1))
          }
          cmatrix(c, cprime) = TLmatrix + Vmatrix;
        }
        
      }else{
        //only include coupling potential
        for(int i = 0; i < N; i++){
          double xi = lagrangeBasis.x(i);
          double ri = a*xi;
          
          for(int j = 0; j < N; j++){
            double xj = lagrangeBasis.x(j);
            double rj = a*xj;
            
            Vmatrix(i,j) = couplingPotential(ri,rj,(*c1),(*c2));
          }
        }
        cmatrix(c, cprime) = Vmatrix;
      }
      cprime++;
    }
    c++;
  }
  return cmatrix;
}

/*calculates the rmatrix from the C matrix*/
arma::cx_mat System::rmatrixCalc(arma::mat<arma::cx_mat> C){
  arma::cx_mat rmatrix(num_channels, num_channels);
  
  std::vector<Channel *>::iterator c1, c2;
  int c = 0, cprime = 0;
  
  for(c1 = channels.begin(); c1 != channels.end(); c1++){
    
    for(c2 = channels.begin(); c2 != channels.end(); c2++){
      
      arma::cx_mat cinv = arma::inv(C(c,cprime));
      double coeff = hbarc*hbarc / (2*a*sqrt((*c1)->getMu() * (*c2)->getMu()));
      
      double expansion = 0;
      for(int i = 0; i < basis_size; i++){
        for(int j = 0; j < basis_size; j++){
          expansion += lagrangeBasis.phi(i,a)*cinv(i,j)*lagrangeBasis.phi(j,a);
        }
      }
      
      rmatrix(c,cprime) = coeff*expansion;
      
      cprime++;
    }
    c++;
  }
  return rmatrix;
}

arma::cx_mat System::umatrixCalc(arma::cx_mat R){
  arma::cx_mat zimatrix(num_channels,num_channels);
  arma::cx_mat zomatrix(num_channels,num_channels);
  arma::cx_double ovalue, ivalue, opvalue, ipvalue, ovalue2, ivalue2, opvalue2, ipvalue2;
  
  std::vector<Channel *>::iterator c1, c2;
  int c = 0, cprime = 0;
  
  for(c1 = channels.begin(); c1 != channels.end(); c1++){
    //calculate wavenumber, rel velocity, and Sommerfeld parameter
    // assuming open channel:
    double kc = sqrt(2*(*c1)->getMu()*(energy - (*c1)->getE())) / hbarc;
    double vc = hbarc*kc / (*c1)->getMu();
    double etac = targ.getZ()*proj.getZ()
      *E2HC/(hbarc*vc);
      
    double hbarx = hbarc*hbarc/(2*mass_unit*(*c1)->getMu());
    double q = sqrt(energy/hbarx);
    
    gsl_sf_result F1,G1,Fp1,Gp1,F2,G2,Fp2,Gp2;
    double exp_F1, exp_G1, exp_F2, exp_G2;
    gsl_sf_coulomb_wave_FG_e(etac,kc*a,(*c1)->getL(),0,&F1,&Fp1,&G1,&Gp1,&exp_F1,&exp_G1);
    ovalue =(G.val+I*F.val);
    ivalue =(G.val-I*F.val);
    opvalue =q*(Gp.val+I*Fp.val);
    ipvalue =q*(Gp.val-I*Fp.val);
    
    for(c2 = channels.begin(); c2 != channels.end(); c2++){
      double kc2 = sqrt(2*(*c2)->getMu()*(energy - (*c2)->getE())) / hbarc;
      double vc2 = hbarc*kc2 / (*c2)->getMu();
      double etac2 = targ.getZ()*proj.getZ()
        *E2HC/(hbarc*vc2);
        
      double hbarx2 = hbarc*hbarc/(2*mass_unit*(*c2)->getMu());
      double q2 = sqrt(energy/hbarx2);
      gsl_sf_coulomb_wave_FG_e(etac2,kc2*a,(*c2)->getL(),0,&F2,&Fp2,&G2,&Gp2,&exp_F2,&exp_G2); 
      ovalue2 =(G2.val+I*F2.val);
      ivalue2 =(G2.val-I*F2.val);
      opvalue2 =q2*(Gp2.val+I*Fp2.val);
      ipvalue2 =q2*(Gp2.val-I*Fp2.val);
      
      zomatrix(c,cprime) = (-1*kc2*a*rmatrix(c,cprime)*opvalue2)
        / (sqrt(kc2*a));
      zimatrix(c,cprime) = (-1*kc2*a*rmatrix(c,cprime)*ipvalue2)
        / (sqrt(kc2*a));
      if(c == cprime){
        zomatrix += ovalue / (sqrt(kc2*a));
        zimatrix += ivalue / (sqrt(kc2*a));
      }
      cprime++;
    }
    c++;
  }
  return arma::inv(zomatrix)*zimatrix;
}

void System::calculateWaveFunction(std::string outFile){
    
  //calculate C matrix
  
  arma::cx_mat cmatrix = cmatrixCalc();
  
  //calculate R matrix
  
  arma::cx_mat rmatrix = rmatrixCalc(cmatrix);
  
  //calculate U matrix
  
  arma::cx_mat umatrix = umatrixCalc(rmatrix);
  
  //calculate wave function
  waveFunction(cmatrix, umatrix);
}
