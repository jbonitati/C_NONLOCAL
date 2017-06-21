#include <cmath>
#include <gsl/gsl_sf.h>
#include <fstream>
#include "system.h"
#include <cstdio>

const double hbarc = 197.3269718;
const double mass_unit=931.494;
const double E2HC = 0.00729927;
const arma::cx_double I(0.0,1.0);

System::System(const double a_size, double e,
  const double m1, const double m2, const double z1, const double z2,
  const int n, std::vector<Channel*> channels_,
  OpticalPotential op, NonLocalOpticalPotential nlop, int basis_size,
  double step, double max, double coupling): 
  a(a_size), energy(e), proj(m1, z1), targ(m2, z2), 
  num_channels(n), channels(channels_),
  pot(op), nlpot(nlop), basis_size(n), step_size(step), r_max(max),
  lagrangeBasis(basis_size), beta(coupling),
  cmatrix(n,n), invcmatrix(n,n), rmatrix(n,n), umatrix(n,n)
{ 
  //calculate C matrix and inverse C matrix
  cmatrixCalc();
  
  //calculate R matrix
  rmatrixCalc();
  
  //calculate U matrix
  umatrixCalc();
}

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
  if(*c1 == *c2){
    return 0;
  }else{
    return beta * (pot.totalPotential(r1, targ, c1) 
      + nlpot.totalPotential(r1,r2,targ,c1));
  }
}

void System::cmatrixCalc(){
  //the C matrix is made up of a different matrix for each pair of channels
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
  //now create inverted C matrix
  for(int i = 0; i < num_channels; i++)
    for(int j = 0; j < num_channels; j++)
      invcmatrix(i,j) = arma::inv(cmatrix(i,j));
}


/*calculates the rmatrix from the inverse C matrix*/
void System::rmatrixCalc(){
  std::vector<Channel *>::iterator c1, c2;
  int c = 0, cprime = 0;
  
  for(c1 = channels.begin(); c1 != channels.end(); c1++){
    
    for(c2 = channels.begin(); c2 != channels.end(); c2++){
      
      arma::cx_mat cinv = Ci->(c,cprime);
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
}

//calculates and returns the collision matrix from the R matrix
void System::umatrixCalc(){
  arma::cx_mat zimatrix(num_channels,num_channels);
  arma::cx_mat zomatrix(num_channels,num_channels);
  arma::cx_double ovalue, ivalue, opvalue, ipvalue, ovalue2, ivalue2, opvalue2, ipvalue2;
  
  std::vector<Channel *>::iterator c1, c2;
  int c = 0, cprime = 0;
  
  for(c1 = channels.begin(); c1 != channels.end(); c1++){
    double kc = (*c1)->getKc(energy);
    (*c1)->io_coulomb_functions(kc*a, energy, targ, proj, 
      &ivalue, &ovalue, &ipvalue, &opvalue);
    
    for(c2 = channels.begin(); c2 != channels.end(); c2++){
      double kc2 = (*c2)->getKc(energy);
      
      (*c2)->io_coulomb_functions(kc2*a, energy, targ, proj, 
        &ivalue2, &ovalue2, &ipvalue2, &opvalue2);
      
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
  umatrix = (arma::inv(zomatrix)*zimatrix);
}

//takes the C and U matrices and stores the wavefunction values for each
//channel in the specified file
void System::waveFunction(std::ofstream file){
  std::cout << "Calculating partial wave functions..." << std::endl;
  file << std::endl << "Printing partial wave functions" << std::endl;
  int c = 0, cprime = 0;
  std::vector<Channel *>::iterator c1, c2;
  
  file << "r";
  for(int i = 1; i <= num_channels; i++){
    file << "\t\tChannel " << i;
  }
  file << std::endl;
  
  for(double r = 0; r < a; r += step_size){
    file << r;
    for(c1 = channels.begin(); c1 != channels.end(); c1++){
      //compute partial wave function u^int_c(c0)(r)
      //sum contains the coefficient calculated over all open channels
      arma::cx_double sum = 0;
      
      arma::cx_double oval, ival, opval, ipval;
      for(c2 = channels.begin(); c2 != channels.end(); c2++){
        //check if channel is open
        
        double 
        
        double kc = (*c2)->getKc(energy);
        double vc = (*c2)->getVc(energy);
        double mu = (*c2)->getMu();
        (*c2)->io_coulomb_functions(kc*a, energy, targ, proj, 
          &ival, &oval, &ipval, &opval);
        double coeff = hbarc*hbarc*kc/(2*mu*sqrt(vc));
        arma::cx_double outersum = 
          -1.0 *coeff*umatrix(cprime, entrance_channel)*opval;
        if(cprime == entrance_channel)
          outersum += coeff*ipval;
        
        arma::cx_double innersum = 0;
        for(int i = 0; i < basis_size; i++){
          for(int j = 0; j < basis_size; j++){
            innersum += lagrangeBasis.phi(i,r)
              *invcmatrix(c,cprime)(i,j)
              *lagrangeBasis.phi(j,a);
          }
        }
        sum += outersum*innersum;
        
        cprime++;
      }
      file << "\t\t" << sum;
      
      c++;
    }
    file << std::endl;
  }
}
