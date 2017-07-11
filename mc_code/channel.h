#pragma once
#include "particle.h"
#include <armadillo>
#include <vector>

class Channel{
  private:
    double E; //channel energy
    int l; //projectile angular momentum
    double j; //total angular momentum
    double m; //projectile spin
    double mu; //reduced mass
    double b; //bloch constant
    double energy; //total energy
    
    double whittaker(double, double, double);//method to calculate general whittaker function
    
  public:
    Channel(const double ec, const int l_, const double j_,
    const double mu_, double energy_, double a, Particle targ, Particle proj):
      E(ec), l(l_), j(j_), m(j_ - l_), mu(mu_), b(0), energy(energy_){ 
      if(energy < ec){
        double kc = getKc();
        b = 2*kc*a*whittaker_prime(2*kc*a, targ, proj)/whittaker(2*kc*a, targ, proj);
      }
    }
    ~Channel(){ }
    
    double getE()const{return E;}
    int getL()const{return l;}
    double getJ()const{return j;}
    double getSpin()const{return m;}
    double getMu()const{return mu;}
    double getB()const{return b;}
    double getKc()const;
    double getVc()const;
    double getEta( Particle targ, Particle proj)const;
    
    double central_potential(double R)const;
    void io_coulomb_functions(double x, Particle targ, Particle proj,
      arma::cx_double *I, arma::cx_double *O, arma::cx_double *Ip, arma::cx_double *Op)const;
    double whittaker(double x, Particle targ, Particle proj);
    double whittaker_prime(double x, Particle targ, Particle proj);  
    
};
