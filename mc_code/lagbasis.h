#include <armadillo>

class LagBasis{
  private:
    const int size;
    const double a; //channel radius
    double * weights;
    double * xi; //zeros of Pn(2x-1)
    arma::vec phi_a;

    double phi(int i, double r) const;
    
    void gauleg(const double x1, const double x2, int n);
    
  public:
    LagBasis(int n, double a);
    ~LagBasis(){
      delete [] weights;
      delete [] xi;
    }
    
    double x(int i)const{ return xi[i];}
    double w(int i)const{ return weights[i];}
    arma::vec get_phi_a()const { return phi_a;}
    arma::rowvec get_phi_r(double r)const;
};
