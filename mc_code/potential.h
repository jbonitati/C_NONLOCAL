#pragma once
#include "channel.h"

class Potential{
  protected:
    //coefficients that can be used by potentials that inherit this one
    double coeff, r, a;
  
  public:
    Potential(const double coefficient, const double radius, const double a_):
      coeff(coefficient), r(radius), a(a_){ }
      
    Potential(const Potential& p): coeff(p.coeff), r(p.r), a(p.a){ }
    
};

class V_Potential : public Potential{
  public:
    double getValue(double R, Particle p)const;
    V_Potential(const double c, const double r_, const double a_):
      Potential(c,r_,a_){ }
    V_Potential(const Potential &p):Potential(p){ }
};

class D_Potential : public Potential{
  public:
    double getValue(double R, Particle p)const;
    D_Potential(const double c, const double r_, const double a_):
      Potential(c,r_,a_){ }
    D_Potential(const Potential &p):Potential(p){ }
};

class SO_Potential : public Potential{
  public:
    double getValue(double R, Particle p, Channel * c)const;
    SO_Potential(const double c, const double r_, const double a_):
      Potential(c,r_,a_){ }
    SO_Potential(const Potential &p):Potential(p){ }
};

class Gauss_potential : public Potential{
  public:
    Gauss_potential(const double coeff, const double exponent): 
      Potential(coeff, exponent, 0){ }
    double getValue(double R) const;
};
  
class CouplingPotential : public Potential{
  private:
    double beta_c;//coupling constant
    double beta_nl; //range of nonlocality
  public:
    CouplingPotential():Potential(0,0,0), beta_c(0), beta_nl(0){ }
    CouplingPotential(double v, double r_, double a_, double bc, double bnl):
      Potential(v,r_,a_), beta_c(bc), beta_nl(bnl){ }
    double getValue(double r1, double r2, Particle p)const;
};
