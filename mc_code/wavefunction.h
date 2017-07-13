#include <vector>
#include <iostream>

class WaveFunction{
  private:
    double r;
    std::vector<double> values;
  public:
    WaveFunction(double r_): r(r_){ }
    void add(double v){ values.push_back(v);}
    
    //function for writing wavefunction values to a stream (as csv)
    friend std::ostream& operator<<(std::ostream &os, const WaveFunction &wf){
      os << wf.r;
      for(std::vector<double>::const_iterator it = wf.values.begin(); it != wf.values.end(); it++)
        os << ", " << (*it);
      os << std::endl;
      return os;
    }
  
};
