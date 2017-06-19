
class Channel{
  private:
    const double E; //energy
    const double l; //projectile angular momentum
    const double j; //total angular momentum
    const double m; //projectile spin
    //spin?
    //mass?
    
  public:
    Channel(const double energy, const double l_, const double j_):
      E(energy), l(l_), j(j_){  }
    
    const double getE()const{return E;}
    const double getL()const{return l;}
    const double getJ()const{return j;}
    const double getSpin()const{return m;}
};
