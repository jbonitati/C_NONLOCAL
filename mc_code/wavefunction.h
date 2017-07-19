#include <vector>
#include <iostream>
#include <armadillo>

using namespace arma;

class WaveFunction{
  private:
    double r;
    std::vector<cx_double> values;
  public:
    WaveFunction(){ }
    WaveFunction(double r_): r(r_){ }
    WaveFunction(double r_, unsigned int num_c): r(r_){ values.reserve(num_c);}
    void add(cx_double v){ values.push_back(v);}
    double getR()const{return r;}
    
    //function for writing wavefunction values to a stream (as csv)
    friend std::ostream& operator<<(std::ostream &os, const WaveFunction &wf){
      os << wf.r;
      for(std::vector<cx_double>::const_iterator it = wf.values.begin(); it != wf.values.end(); it++)
        os << ", " << std::abs((*it));
      os << std::endl;
      return os;
    }
  
};
