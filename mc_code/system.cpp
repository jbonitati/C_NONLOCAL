#include <cmath>
#include <iostream>
#include <vector>
#include <boost/timer/timer.hpp>
#include <boost/math/special_functions/round.hpp>

#include "system.h"

using namespace boost::numeric::ublas; //matrix
using boost::math::iround; //iround
using namespace arma; //cx_double

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
  //std::cout << cmatrix;
  
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
  
  /*
  std::cout << "The following matrix should be identity if using real potential:" << std::endl;
  std::cout << umatrix.t()*umatrix << std::endl;
  std::cout << "u_int(a):";
  for(unsigned int i = 0; i < c.size(); i++) std::cout << internalWaveFunction(i,a);
  std::cout << "\nu_ext(a):";
  for(unsigned int i = 0; i < c.size(); i++) std::cout << externalWaveFunction(i,a);
  std::cout << std::endl;
  double h = 0.01;
  std::cout << "u_int'(a):";
  for(unsigned int i = 0; i < c.size(); i++){
    std::cout << 
    (internalWaveFunction(i,a-h) + internalWaveFunction(i,a))/h;
  }
  std::cout << "\nu_ext'(a):";
  for(unsigned int i = 0; i < c.size(); i++){
    std::cout << 
    (externalWaveFunction(i,a-h) + externalWaveFunction(i,a))/h;
  }
  std::cout << std::endl;
  */
}

//returns the coupling potential between two channels
double System::couplingPotential(double r1, double r2,
  unsigned int c1, unsigned int c2){
  //return 0;
  return coupling_matrix(c1,c2).getValue(r1,r2,targ);
}

void System::cmatrixCalc()
{
  boost::timer::auto_cpu_timer t("C matrix time: %w sec real, %t sec CPU\n");
  //the C matrix is made up of a different matrix for each pair of channels
  //i.e. it is (c*N)x(c*N) where index (c*N+i,c'*N+j) corresponds to 
  // basis function i in channel c and basis function j in channel c'
  int N = basis_size;
  
  //TLmatrix uses the T and L operators
  //Vmatrix uses the potential, centrifugal, and energy operators
  for(unsigned int c = 0; c < channels.size(); c++)
  {
    Channel * c1 = &channels[c];
    
    for(unsigned int cprime = 0; cprime < channels.size(); cprime++)
    {
      //Channel * c2 = &channels[cprime];
      arma::cx_mat Vmatrix(N,N, arma::fill::zeros);
      arma::cx_mat TLmatrix(N, N, arma::fill::zeros);

      if(c == cprime)
      {
        //include TLmatrix as well as Vmatrix for c = c'
    
        double mu = c1->getMu();
        double TLcoeff = hbarc*hbarc/(2*mass_unit*mu);
        
        for(int i = 0; i < N; i++)
        { 
          double xi = lagrangeBasis.x(i);
          double ri = a*xi;
          
          TLmatrix(i,i) = TLcoeff*((4*N*N+4*N+3)*xi*(1-xi)-6*xi+1)/
            (3*a*a*xi*xi*(1-xi)*(1-xi)) - c1->getB();

          Vmatrix(i,i) = c1->central_potential(ri) 
            + c1->getE() - energy
            + pot.totalPotential(ri, targ, c1)
            + proj.coulomb_potential_to(ri, coulomb_radius, targ);
          for(int j = 0; j < N; j++)
          {
            double xj = lagrangeBasis.x(j);
            double rj = a*xj;
            
            if(j > i)
            {
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
        
      }
      else
      {
        //only include coupling potential for different channels
        for(int i = 0; i < N; i++)
        {
          double xi = lagrangeBasis.x(i);
          double ri = a*xi;
          
          for(int j = 0; j < N; j++)
          {
            double xj = lagrangeBasis.x(j);
            double rj = a*xj;
            
            Vmatrix(i,j) = a*sqrt(lagrangeBasis.w(i)*lagrangeBasis.w(j))*
              couplingPotential(ri,rj,c,cprime);
          }
        }
      }

      cmatrix.submat(c*N, cprime*N, (c+1)*N-1, (cprime+1)*N-1) = TLmatrix + Vmatrix;
    }
  }
}

/*calculates the inverse of the C matrix*/
void System::invertCmatrix()
{
  boost::timer::auto_cpu_timer t("C matrix inversion time: %w sec real, %t sec CPU\n");
  invcmatrix = arma::inv(cmatrix);
}

/*calculates the rmatrix from the inverse C matrix*/
void System::rmatrixCalc()
{
  boost::timer::auto_cpu_timer t("R matrix time: %w sec real, %t sec CPU\n");
  
  for(unsigned int c = 0; c < channels.size(); c++)
  {
    Channel * c1 = &channels[c];
    for(unsigned int cprime = 0; cprime < channels.size(); cprime++)
    {
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
void System::umatrixCalc()
{
  boost::timer::auto_cpu_timer t("U matrix time: %w sec real, %t sec CPU\n");
  arma::cx_mat zimatrix(channels.size(),channels.size(), arma::fill::zeros);
  arma::cx_mat zomatrix(channels.size(),channels.size(), arma::fill::zeros);
  cx_double ovalue, ivalue, opvalue, ipvalue;//, ovalue2, ivalue2, opvalue2, ipvalue2;
  
  for(unsigned int c = 0; c < channels.size(); c++)
  {
    //Channel * c1 = &channels[c];
    //double kc = c1->getKc();
    //c1->io_coulomb_functions(kc*a, energy, targ, proj, 
    //  &ivalue, &ovalue, &ipvalue, &opvalue);
    
    for(unsigned int cprime = 0; cprime < channels.size(); cprime++)
    {
      Channel * c2 = &channels[cprime];
      double kc2 = c2->getKc();
      cx_double coeff = 1.0 / (sqrt(kc2*a)); 
      
      c2->io_coulomb_functions(kc2*a, targ, proj, 
        &ivalue, &ovalue, &ipvalue, &opvalue);
      
      //old python program did not have k in the following equations
      //this is due to the way it took the derivative differently
      zomatrix(c,cprime) = coeff*(-1*kc2*
        a*rmatrix(c,cprime)*opvalue);
      zimatrix(c,cprime) = coeff*(-1*kc2*
        a*rmatrix(c,cprime)*ipvalue);
        
      if(c == cprime)
      {
        zomatrix(c,cprime) += coeff*ovalue;
        zimatrix(c,cprime) += coeff*ivalue;
      }
    }
  }
  //std::cout << zomatrix << zimatrix;
  umatrix = (arma::inv(zomatrix)*zimatrix);
  for(unsigned int c = 0; c < channels.size(); c++)
  {
    for(unsigned int cprime = 0; cprime < channels.size(); cprime++)
      umatrix(c,cprime) *= sqrt(channels[cprime].getKc() / channels[c].getKc());
  }
}

//Calculates the wave function for channel c at value r
cx_double System::internalWaveFunction(unsigned int c, double r)
{

  //Channel * c1 = &channels[c];
  //double vc1 = c1->getVc();
  //double kc0 = channels[entrance_channel].getKc();
  
  //compute partial wave function u^int_c(c0)(r)
  cx_double wfvalue = 0;
  cx_double oval, ival, opval, ipval;
  
  for(unsigned int cprime = 0;
    cprime < channels.size();
    cprime++)
  {
    if( channels[cprime].isOpen())
    {
    Channel * c2 = &channels[cprime];
  
    double kc = c2->getKc();
    //double vc = c2->getVc();
    double mu = c2->getMu();
    c2->io_coulomb_functions(kc*a, targ, proj, 
      &ival, &oval, &ipval, &opval);
    
    cx_double coeff = nrmlz*hbarc*hbarc*kc
      /(2*mass_unit*mu);//*sqrt(channels[c].getVc()));
    cx_double outersum = 
      -1.0 *coeff*umatrix(cprime, entrance_channel)*opval;
    if(cprime == entrance_channel)
      outersum += coeff*ipval;
    //if(c2->getB() != 0)
    //  outersum += -1.0*hbarc*hbarc/(2*mass_unit*mu)*
    //    c2->getB()/a * externalWaveFunction(cprime,a);

    arma::rowvec phir = lagrangeBasis.get_phi_r(r);
    arma::cx_mat cinv = invcmatrix.submat(c*basis_size, cprime*basis_size, 
        (c+1)*basis_size - 1, (cprime + 1)*basis_size - 1);
    arma::vec phia = lagrangeBasis.get_phi_a();
    cx_double innersum = arma::as_scalar(phir*cinv*phia);
    
    wfvalue += outersum*innersum;
    }
  }
  return wfvalue;
}

//Calculates the wave function for channel c at position r
cx_double System::externalWaveFunction(unsigned int c, double r)
{
  //double k0 = channels[entrance_channel].getKc();
  WaveFunction wf(r, channels.size());
  cx_double oval, ival, opval, ipval;
  Channel * c1 = &channels[c];
  double kc = c1->getKc();
  //double vc = c1->getVc();
  //double etac = c1->getEta(energy, targ, proj);
  //double mu = c1->getMu();
  c1->io_coulomb_functions(kc*r, targ, proj, 
    &ival, &oval, &ipval, &opval);
  cx_double wfvalue;
  if(c1->isOpen())
  {
    //Open Channel
    double coeff = 1.0;///sqrt(c1->getVc());
    wfvalue = -1.0 *coeff*umatrix(c, entrance_channel)*oval;
    if(c == entrance_channel)
      wfvalue += coeff*ival;
      
  }
  else
  {
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
  return wfvalue;
}

//calculates the wavefunction for all channels at distance r
//Returns the wave function in WaveFunction object, which
// contains the distance  r  and values for each channel
WaveFunction System::calculateWaveFunction(double r)
{
  WaveFunction wf(r);
  for(unsigned int c = 0; c < channels.size(); c++)
  {
    if(r < a)
    {
      wf.add(internalWaveFunction(c,r));
    }
    else
    {
      wf.add(externalWaveFunction(c,r));
    }
  }
  return wf;
}

//calculates the wave functions for each channel and stores them in file
void System::waveFunctions(boost::filesystem::ofstream& file)
{
  boost::timer::auto_cpu_timer t("Wave Function calculation time: %w sec real, %t sec CPU\n");
  int num_values = iround(r_max / step_size) + 1;
  std::vector<WaveFunction> wfvalues(num_values);
  
  int i;
  #pragma omp parallel for
  for(i = 0; i < num_values; i++)
  {
    wfvalues[i] = calculateWaveFunction(i*step_size);
  }
  
  std::cout << "Printing Wave functions to file..." << std::endl;
  
  //file << "r";
  //for(unsigned int i = 1; i <= channels.size(); i++)
  //{
    //file << " Channel_" << i << "(Re)";
    //file << " Channel_" << i << "(Im)";
  //  file << " Channel_" << i << "(abs)";
  //}
  //file << std::endl;
  
  for(std::vector<WaveFunction>::iterator it = wfvalues.begin(); it!= wfvalues.end(); it++)
  {
    file << (*it);
  }  
}

//calculates the total potential at all points used in wave function calculations
// and stores them in file
void System::plotPotential(boost::filesystem::ofstream& file)
{
  file << "r Re(potential) Im(potential) Coupling" <<std::endl;
  for(double r = step_size; r <= r_max; r+= step_size)
  {
    file << r << " " 
      << std::real(pot.totalPotential(r, targ, &channels[0]) 
        + nlpot.totalPotential(r,r,targ,&channels[0]))
      << " "
      << std::imag(pot.totalPotential(r, targ, &channels[0]) 
        + nlpot.totalPotential(r,r,targ,&channels[0]))
      << " "
      << couplingPotential(r,r,0,channels.size()-1)
      << std::endl;
  }
  
}
