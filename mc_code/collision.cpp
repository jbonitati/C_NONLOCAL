
class Collision{
  private:
    const double a; // channel radius
    double b; //Bloch constant
    
    const Particle proj;  //projectile particle
    const Particle targ; //target particle
    
    Channel* channels;
    
    
    
  public:
    Collision~(){
      delete [] channels;
    }
    
    Collision(const double a_size, double bloch, 
      const double m1, const double m2, const double z1, const double z2,
      Channel * channels_): 
      a(a_size), b(bloch), proj(m1, z1), targ(m2, z2), channels(channels_)
    { }
    
    const double getA() const{return a;}
    double getB() const{return b;}
    const Particle getProj() const{return proj;}
    const Particle getTarg() const{return targ;}
    Channel* getChannels() const{return channels;}
}
