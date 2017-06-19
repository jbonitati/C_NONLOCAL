
#include <cmath>
const std::complex<double> I(0., 1.);

std::complex<double> OpticalPotential::Volume_potential
(double radius, Particle p){
  return Volume_Re.getValue(radius, p) 
    + I*Volume_Im.getValue(radius, p);
}

std::complex<double> OpticalPotential::Surface_potential
(double radius, Particle p){
  return Surface_Re.getValue(radius, p) 
    + I*Surface_Im.getValue(radius, p);
}

std::complex<double> OpticalPotential::Spin_orbit_potential
(double radius, Particle p, Channel * c){
  return Surface_Re.getValue(radius, p, c)
    + I*Surface_Im.getValue(radius, p, c);
}

virtual std::complex<double> OpticalPotential::totalPotential
(double radius, Particle p, Channel * c){
  return Volume_potential(radius, p) 
    + Surface_potential(radius, p) 
    + Spin_orbit_potential(radius, p, c);
}

std::complex<double> NonLocalOpticalPotential::totalPotential
(double r1, double r2, Particle p, Channel  * c){
  double ravg = (r1+r2)/2.0;
  
  std::complex<double> local_part 
    = OpticalPotential::totalPotential(ravg, p, c);
    
  std::complex<double> nonlocal_part 
    = 1.0/(pow(beta,3.0)*pow(PI,3.0/2.0))*exp(-pow((r1-r2)/beta,2.0));
    
  return local_part * nonlocal_part;
}
