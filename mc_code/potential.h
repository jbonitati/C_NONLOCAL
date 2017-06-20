#include "channel.h"
#include <stdexcept>

class Potential{
  private:
    const double coeff, r, a;
  
  public:
    Potential(const double coefficient, const double radius, const double a_):
      coeff(coefficient), r(radius), a(a_){ }
      
    Potential(const Potential& p){
      coeff = p.coeff;
      r = p.r;
      a = p.a;
    }
    
};

class V_Potential : public Potential{
  public:
    double getValue(double R, Particle p);
    V_Potential(const double c, const double r_, const double a_):
      Potential(c,r_,a_){ }
    V_Potential(const V_Potential &p):Potential(p){ }
};

class D_Potential : public Potential{
  public:
    double getValue(double R, Particle p);
    D_Potential(const double c, const double r_, const double a_):
      Potential(c,r_,a_){ }
    D_Potential(const D_Potential &p):Potential(p){ }
};

class SO_Potential : public Potential{
  public:
    double getValue(double R, Particle p, Channel * c);
    SO_Potential(const double c, const double r_, const double a_):
      Potential(c,r_,a_){ }
    SO_Potential(const SO_Potential &p):Potential(p){ }
};
