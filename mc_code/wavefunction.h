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
    
    //function for writing wavefunction values to a stream
    friend std::ostream& operator<<(std::ostream &os, const WaveFunction &wf);
};
