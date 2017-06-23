#include "particle.h"
#include <armadillo>
#pragma once

class Channel{
  private:
    double E; //energy
    int l; //projectile angular momentum
    double j; //total angular momentum
    double m; //projectile spin
    double mu; //reduced mass
    double b; //bloch constant
    
  public:
    Channel(const double energy, const int l_, const double j_,
    const double mu_):
      E(energy), l(l_), j(j_), m(j_ - l_), mu(mu_), b(0){ }
    
    double getE()const{return E;}
    int getL()const{return l;}
    double getJ()const{return j;}
    double getSpin()const{return m;}
    double getMu()const{return mu;}
    double getB()const{return b;}
    double getKc(double energy)const;
    double getVc(double energy)const;
    double getEta(double energy, Particle targ, Particle proj)const;
    
    double central_potential(double R)const;
    void io_coulomb_functions(double x, double energy, Particle targ, Particle proj,
      arma::cx_double *I, arma::cx_double *O, arma::cx_double *Ip, arma::cx_double *Op)const;
};
