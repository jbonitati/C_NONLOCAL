#include "particle.h"
#include <armadillo>
#pragma once

class Channel{
  private:
    const double E; //energy
    const double l; //projectile angular momentum
    const double j; //total angular momentum
    const double m; //projectile spin
    double mu; //reduced mass
    double b; //bloch constant
    //spin?
    //mass?
    
  public:
    Channel(const double energy, const double l_, const double j_,
    const double m_):
      E(energy), l(l_), j(j_), m(m_), b(0){  }
    
    const double getE()const{return E;}
    const double getL()const{return l;}
    const double getJ()const{return j;}
    const double getSpin()const{return m;}
    double getMu()const{return mu;}
    const double getB()const{return b;}
    double getKc(double energy);
    double getVc(double energy);
    
    double central_potential(double R);
    void io_coulomb_functions(double x, double energy, Particle targ, Particle proj,
      arma::cx_double *I, arma::cx_double *O, arma::cx_double *Ip, arma::cx_double *Op);
};
