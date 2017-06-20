
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
    const double m_)):
      E(energy), l(l_), j(j_), m(m_), b(0){  }
    
    const double getE()const{return E;}
    const double getL()const{return l;}
    const double getJ()const{return j;}
    const double getSpin()const{return m;}
    double getMu()const{return mu;}
    const double getB()const{return b;}
    
    double central_potential(double R, Particle proj, Particle targ);
    double coulomb_potential(double R, Particle proj, Particle targ);
};
