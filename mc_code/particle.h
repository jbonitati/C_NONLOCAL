#pragma once

class Particle{
  private:
    const double m; //mass number
    const double z; //proton number
    const double n; //neutron number
  
  public:
    Particle(const double mass, const double protons): 
      m(mass), z(protons), n(mass-protons){ }
    const double getM() const {return m;}
    const double getZ() const {return z;}
    const double getN() const {return n;}
    
    double coulomb_potential_to(double R, Particle targ)const;
};
