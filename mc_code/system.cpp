#include <cmath>
#include <gsl/gsl_sf.h>
#include <fstream>
#include "system.h"
#include <iostream>
#include <vector>
#include "wavefunction.h"
#include <complex>
#include "potential.h"
#include "constants.h"

using namespace boost::numeric::ublas;
using namespace arma; //cx_double


System::System(const double a_size, double e,
  Particle projectile, Particle target,
  int c0, std::vector<Channel> c,
  OpticalPotential op, NonLocalOpticalPotential nlop, int b_size,
  double step, double max, matrix<CouplingPotential> coupling): 
  a(a_size), energy(e), proj(projectile), targ(target), 
  entrance_channel(c0), channels(c),
  pot(op), nlpot(nlop), basis_size(b_size), step_size(step), r_max(max),
  lagrangeBasis(basis_size, a_size), cmatrix(c.size()*basis_size, c.size()*basis_size),
  invcmatrix(c.size()*basis_size, c.size()*basis_size), rmatrix(c.size(), c.size()), 
  umatrix(c.size(),c.size()), coupling_matrix(coupling)
{
  
  std::cout << "Calculating C matrix..." << std::endl;
  //calculate C matrix and inverse C matrix
  cmatrixCalc();
  
  std::cout << "Calculating R matrix..." << std::endl;
  //calculate R matrix
  rmatrixCalc();
  std::cout << rmatrix << std::endl;
  
  std::cout << "Calculating Collision (U) matrix..." << std::endl;
  //calculate U matrix
  umatrixCalc();
  std::cout << umatrix << std::endl;
  
  //std::cout << channels[0].getKc() << " " << a_size << std::endl;
  //channels[0].coulomb_test(channels[0].getKc()*a_size, targ, proj);
}


//calculates the coupling potential between two channels
double System::couplingPotential(double r1, double r2, int c1, int c2){
  return coupling_matrix(c1,c2).getValue(r1,r2,targ);
}

void System::cmatrixCalc(){
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
            + proj.coulomb_potential_to(ri, targ);
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
            
            Vmatrix(i,j) = a*sqrt(lagrangeBasis.w(i)*lagrangeBasis.w(j))
              *couplingPotential(ri,rj,c,cprime);
          }
        }
        cmatrix.submat(c*N, cprime*N, (c+1)*N-1, (cprime+1)*N-1) = Vmatrix;
      }
    }
  }
  std::cout << "Calculating inverse C matrix..." << std::endl;
  //std::cout << arma::inv(arma::real(cmatrix(0,1))) << std::endl;
  //now create inverted C matrix
  //for(unsigned int i = 0; i < channels.size(); i++){
  //  for(unsigned int j = 0; j < channels.size(); j++){
  invcmatrix = arma::inv(cmatrix);
  //  }
  //}
}


/*calculates the rmatrix from the inverse C matrix*/
void System::rmatrixCalc(){
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
  arma::cx_mat zimatrix(channels.size(),channels.size());
  arma::cx_mat zomatrix(channels.size(),channels.size());
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
        zomatrix += coeff*ovalue;
        zimatrix += coeff*ivalue;
      }
    }
  }
  umatrix = (arma::inv(zomatrix)*zimatrix);
}

//takes the C and U matrices and stores the wavefunction values for each
//channel in the specified file
void System::waveFunction(boost::filesystem::ofstream& file){
  
  std::vector<WaveFunction> wfvalues;
  wfvalues.reserve(std::floor(r_max / step_size) + 1);
  
  double r;
  std::cout << 0;
  for(r = 0; r < a && r <= r_max; r += step_size){ 
    std::cout << "\rradius: " << r << std::flush;   
    //file << r;
    WaveFunction wf(r);
    for(unsigned int c = 0; c < channels.size(); c++){
      Channel * c1 = &channels[c];
      double vc1 = c1->getVc();
      
      //compute partial wave function u^int_c(c0)(r)
      //sum contains the coefficient calculated over all open channels
      cx_double sum = 0;
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
          
        //if(r == 0){
          /*std::cout << "I(ka) " << ival << ", O(ka) " << oval << ", I'(ka) " << ipval << ", O'(ka) " << opval << std::endl;
          std::cout << kc << std::endl;
          std::cout << umatrix(cprime, entrance_channel) << std::endl;
          std::cout << nrmlz*(ipval - umatrix(cprime, entrance_channel)*opval) << std::endl;
          */
          //std::cout << nrmlz*outersum << std::endl;
        //}
        
        cx_double innersum = 0;
        arma::rowvec phir = lagrangeBasis.get_phi_r(r);
        arma::cx_mat cinv = invcmatrix.submat(c*basis_size, cprime*basis_size, 
            (c+1)*basis_size - 1, (cprime + 1)*basis_size - 1);
        arma::vec phia = lagrangeBasis.get_phi_a();
        innersum = arma::as_scalar(phir*cinv*phia);
        /*for(int i = 0; i < basis_size; i++){
          for(int j = 0; j < basis_size; j++){
            innersum += lagrangeBasis.phi(i+1,r)
              *invcmatrix(c*basis_size + i, cprime*basis_size + j)
              *lagrangeBasis.phi(j+1,a);
          }
        }*/
        sum += outersum*innersum;
        
      }
      //sum *= nrmlz;
      //file << ", " << std::real(sum);
      wf.add(std::real(sum));
    }
    wfvalues.push_back(wf);
  }
  
  //now print the external wave function calculated from coulomb scattering
  //with the collision matrix
  while(r <= r_max){
    std::cout << "\rradius: " << r << std::flush;
   //file << r;
   WaveFunction wf(r);
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
        
        double coeff = 1.0/sqrt(vc);
        wfvalue = -1.0 *coeff*umatrix(c, entrance_channel)*oval;
        if(c == entrance_channel)
          wfvalue += coeff*ival;
          
      }else{
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
      //file << ", " << std::real(wfvalue);
      wf.add(std::real(wfvalue));
    }
    wfvalues.push_back(wf);
  
   r += step_size; 
  }
  std::cout << std::endl;
  std::cout << "Printing Wave functions to file..." << std::endl;
  
  file << "r";
  //file << 
  for(unsigned int i = 1; i <= channels.size(); i++){
    file << ", Channel " << i;
  }
  file << std::endl;
  for(std::vector<WaveFunction>::iterator it = wfvalues.begin(); it != wfvalues.end(); it++){
    file << (*it);
    //file << (*it).getR() << ", " << pot.totalPotential((*it).getR(), targ, &channels[0]) << std::endl;
  }
}
