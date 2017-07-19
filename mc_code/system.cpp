#include <cmath>
#include <gsl/gsl_sf.h>
#include "system.h"
#include <iostream>
#include <vector>
#include <boost/timer/timer.hpp>
#include <boost/math/special_functions/round.hpp>

using namespace boost::numeric::ublas;
using boost::math::iround;
using namespace arma; //for cx_double typedef


System::System(const double a_size, double e,
  Particle projectile, Particle target, double rc,
  int c0, std::vector<Channel> c,
  OpticalPotential op, NonLocalOpticalPotential nlop, int b_size,
  double step, double max, matrix<CouplingPotential> coupling): 
  a(a_size), energy(e), proj(projectile), targ(target), coulomb_radius(rc),
  entrance_channel(c0), channels(c),
  pot(op), nlpot(nlop), basis_size(b_size), step_size(step), r_max(max),
  lagrangeBasis(basis_size, a_size), cmatrix(c.size()*basis_size, c.size()*basis_size),
  invcmatrix(c.size()*basis_size, c.size()*basis_size), rmatrix(c.size(), c.size()), 
  umatrix(c.size(),c.size()), coupling_matrix(coupling)
{
  
  std::cout << "Calculating C matrix..." << std::endl;
  cmatrixCalc();
  
  std::cout << "Calculating inverse C matrix..." << std::endl;
  invertCmatrix();
  
  std::cout << "Calculating R matrix..." << std::endl;
  rmatrixCalc();
  //if(channels.size() <= 4)
    std::cout << rmatrix << std::endl;
  
  std::cout << "Calculating Collision (U) matrix..." << std::endl;
  umatrixCalc();
  //if(channels.size() <= 4)
    std::cout << umatrix << std::endl;
  
  //std::cout << "The following matrix should be identity if using real potential:" << std::endl;
  //std::cout << umatrix.t()*umatrix << std::endl;
  
}


//returns the coupling potential between two channels
double System::couplingPotential(double r1, double r2,
  unsigned int c1, unsigned int c2){
  //return 0;
  return coupling_matrix(c1,c2).getValue(r1,r2,targ);
}

void System::cmatrixCalc(){
  boost::timer::auto_cpu_timer t;
  //the C matrix is made up of a different matrix for each pair of channels
  int N = basis_size;
  
  //TLmatrix uses the T and L operators
  //Vmatrix uses the potential, centrifugal, and energy operators
  for(unsigned int c = 0; c < channels.size(); c++){
    Channel * c1 = &channels[c];
    
    for(unsigned int cprime = 0; cprime < channels.size(); cprime++){
      //Channel * c2 = &channels[cprime];
      arma::cx_mat Vmatrix(N,N, arma::fill::zeros);
      if(c == cprime){
        //include TLmatrix as well as Vmatrix for c = c'
        arma::cx_mat TLmatrix(N, N, arma::fill::zeros);
    
        double mu = c1->getMu();
        double TLcoeff = hbarc*hbarc/(2*mass_unit*mu);
        
        for(int i = 0; i < N; i++){ 
          double xi = lagrangeBasis.x(i);
          double ri = a*xi;
          
          TLmatrix(i,i) = TLcoeff*((4*N*N+4*N+3)*xi*(1-xi)-6*xi+1)/
            (3*a*a*xi*xi*(1-xi)*(1-xi)) - c1->getB();
          Vmatrix(i,i) = c1->central_potential(ri) 
            + c1->getE() - energy
            + pot.totalPotential(ri, targ, c1)
            + proj.coulomb_potential_to(ri, coulomb_radius, targ);
          for(int j = 0; j < N; j++){
            double xj = lagrangeBasis.x(j);
            double rj = a*xj;
            
            if(j > i){
              TLmatrix(i,j)=TLcoeff 
              * pow(-1,i+j)/(a*a*sqrt(xi*xj*(1-xi)*(1-xj)))
              * (N*N + N + 1 + (xi+xj-2*xi*xj) / ((xi-xj)*(xi-xj)) 
                - 1.0/(1-xi) - 1.0/(1-xj)) - c1->getB();
              TLmatrix(j,i)=TLmatrix(i,j);
            }
            
            Vmatrix(i,j) += a*sqrt(lagrangeBasis.w(i)*lagrangeBasis.w(j)) 
              * nlpot.totalPotential(ri, rj, targ, c1);
          }
        }
        cmatrix.submat(c*N, cprime*N, (c+1)*N-1, (cprime+1)*N-1) = TLmatrix + Vmatrix;
        
      }else{
        //only include coupling potential for different channels
        for(int i = 0; i < N; i++){
          double xi = lagrangeBasis.x(i);
          double ri = a*xi;
          
          for(int j = 0; j < N; j++){
            double xj = lagrangeBasis.x(j);
            double rj = a*xj;
            
            Vmatrix(i,j) = a*sqrt(lagrangeBasis.w(i)*lagrangeBasis.w(j))*
              couplingPotential(ri,rj,c,cprime);
          }
        }
        cmatrix.submat(c*N, cprime*N, (c+1)*N-1, (cprime+1)*N-1) = Vmatrix;
      }
    }
  }
}

/*calculates the inverse of the C matrix*/
void System::invertCmatrix(){
  boost::timer::auto_cpu_timer t;
  invcmatrix = arma::inv(cmatrix);
}

/*calculates the rmatrix from the inverse C matrix*/
void System::rmatrixCalc(){
  boost::timer::auto_cpu_timer t;
  for(unsigned int c = 0; c < channels.size(); c++){
    Channel * c1 = &channels[c];
    for(unsigned int cprime = 0; cprime < channels.size(); cprime++){
      Channel * c2 = &channels[cprime];
      arma::cx_mat cinv = invcmatrix.submat(c*basis_size,cprime*basis_size,
        (c+1)*basis_size-1, (cprime+1)*basis_size-1);
      
      double coeff = hbarc*hbarc / (2*a*mass_unit*sqrt(c1->getMu() * c2->getMu()));
      arma::vec phia = lagrangeBasis.get_phi_a();
      cx_double expansion = arma::as_scalar(phia.t()*cinv*phia);
    
      rmatrix(c,cprime) = coeff*expansion;
      
    }
  }
}

//calculates and returns the collision matrix from the R matrix
void System::umatrixCalc(){
  boost::timer::auto_cpu_timer t;
  arma::cx_mat zimatrix(channels.size(),channels.size(), arma::fill::zeros);
  arma::cx_mat zomatrix(channels.size(),channels.size(), arma::fill::zeros);
  cx_double ovalue, ivalue, opvalue, ipvalue;//, ovalue2, ivalue2, opvalue2, ipvalue2;
  
  for(unsigned int c = 0; c < channels.size(); c++){
    //Channel * c1 = &channels[c];
    //double kc = c1->getKc();
    //c1->io_coulomb_functions(kc*a, energy, targ, proj, 
    //  &ivalue, &ovalue, &ipvalue, &opvalue);
    
    for(unsigned int cprime = 0; cprime < channels.size(); cprime++){
      Channel * c2 = &channels[cprime];
      double kc2 = c2->getKc();
      cx_double coeff = 1.0 / (sqrt(kc2*a)); 
      
      c2->io_coulomb_functions(kc2*a, targ, proj, 
        &ivalue, &ovalue, &ipvalue, &opvalue);
      
      //python program did not have k in the following equations
      //this is due to the way it took the derivative differently
      zomatrix(c,cprime) = coeff*(-1*kc2*
        a*rmatrix(c,cprime)*opvalue);
      zimatrix(c,cprime) = coeff*(-1*kc2*
        a*rmatrix(c,cprime)*ipvalue);
        
      if(c == cprime){
        zomatrix(c,cprime) += coeff*ovalue;
        zimatrix(c,cprime) += coeff*ivalue;
      }
    }
  }
  //std::cout << zomatrix << zimatrix;
  umatrix = (arma::inv(zomatrix)*zimatrix);
}

//Calculates the wave function for all channels at r < a
WaveFunction System::internalWaveFunction(double r){
  WaveFunction wf(r, channels.size());
  for(unsigned int c = 0; c < channels.size(); c++){
    Channel * c1 = &channels[c];
    double vc1 = c1->getVc();
    
    //compute partial wave function u^int_c(c0)(r)
    //sum contains the coefficient calculated over all open channels
    cx_double wfvalue = 0;
    cx_double oval, ival, opval, ipval;
    
    for(unsigned int cprime = 0; 
      cprime < channels.size() && channels[cprime].isOpen();
      cprime++)
    {
      Channel * c2 = &channels[cprime];
    
      double kc = c2->getKc();
      //double vc = c2->getVc();
      double mu = c2->getMu();
      c2->io_coulomb_functions(kc*a, targ, proj, 
        &ival, &oval, &ipval, &opval);
      
      cx_double coeff = nrmlz*hbarc*hbarc*kc
        /(2*mass_unit*mu*sqrt(vc1));
      cx_double outersum = 
        -1.0 *coeff*umatrix(cprime, entrance_channel)*opval;
      if(cprime == entrance_channel)
        outersum += coeff*ipval;

      arma::rowvec phir = lagrangeBasis.get_phi_r(r);
      arma::cx_mat cinv = invcmatrix.submat(c*basis_size, cprime*basis_size, 
          (c+1)*basis_size - 1, (cprime + 1)*basis_size - 1);
      arma::vec phia = lagrangeBasis.get_phi_a();
      cx_double innersum = arma::as_scalar(phir*cinv*phia);
      wfvalue += outersum*innersum;
      
    }
    wf.add(wfvalue);
  }
  return wf;
}

//Calculates the wave function for all channels at r > a
WaveFunction System::externalWaveFunction(double r){
  WaveFunction wf(r, channels.size());
  for(unsigned int c = 0; c < channels.size(); c++){
    cx_double oval, ival, opval, ipval;
    Channel * c1 = &channels[c];
    double kc = c1->getKc();
    double vc = c1->getVc();
    //double etac = c1->getEta(energy, targ, proj);
    //double mu = c1->getMu();
    c1->io_coulomb_functions(kc*r, targ, proj, 
      &ival, &oval, &ipval, &opval);
    cx_double wfvalue;
    if(energy >= c1->getE()){
      //Open Channel
      
      double coeff = 1.0/sqrt(vc);
      wfvalue = -1.0 *coeff*umatrix(c, entrance_channel)*oval;
      if(c == entrance_channel)
        wfvalue += coeff*ival;
        
    }else{
      //Closed Channel
      
      double coeff = 1.0/c1->whittaker(2*kc*a, targ, proj);
      cx_double sum = 0;
      for(unsigned int cprime = 0; 
        cprime < channels.size() && channels[cprime].isOpen();
        cprime++)
      {
        Channel * c2 = &channels[cprime];
        cx_double oval2, ival2, opval2, ipval2;
        double kc2 = c2->getKc();
        double mu2 = c2->getMu();
        c2->io_coulomb_functions(kc2*a, targ, proj, 
          &ival2, &oval2, &ipval2, &opval2);
        cx_double sumcoeff = sqrt(mass_unit*mu2*kc2/hbarc)*a*rmatrix(c,cprime);
        sum += -1.0*sumcoeff*umatrix(cprime,entrance_channel)*opval2;
        if(cprime == entrance_channel)
          sum += sumcoeff*ipval2;
      }
      wfvalue = coeff*sum*c1->whittaker(2*kc*r, targ, proj);
    }
    wfvalue *= nrmlz;
    wf.add(wfvalue);
  }
  return wf;
}

//calculates the wavefunction for all channels at distance r
//Returns the wave function in WaveFunction object, which
// contains the distance  r  and values for each channel
WaveFunction System::calculateWaveFunction(double r){
  if(r < a){
    return internalWaveFunction(r);
  }else{
    return externalWaveFunction(r);
  }
}

//calculates the wave functions for each channel and stores them in file
void System::waveFunctions(boost::filesystem::ofstream& file){
  boost::timer::auto_cpu_timer t;
  int num_values = iround(r_max / step_size) + 1;
  std::vector<WaveFunction> wfvalues(num_values);
  
  int i;
  #pragma omp parallel for
  for(i = 0; i < num_values; i++){
    wfvalues[i] = calculateWaveFunction(i*step_size);
  }
  
  std::cout << std::endl;
  std::cout << "Printing Wave functions to file..." << std::endl;
  
  file << "r";
  for(unsigned int i = 1; i <= channels.size(); i++){
    file << ", Channel " << i;
  }
  file << std::endl;
  for(std::vector<WaveFunction>::iterator it = wfvalues.begin(); it!= wfvalues.end(); it++){
    file << (*it);
  }
}
