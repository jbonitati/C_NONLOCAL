
class LagBasis{
  private:
    const int size;
    const double a; //channel radius
    double * weights;
    double * xi; //zeros of Pn(2x-1)
    
    void gauleg(const double x1, const double x2, int n);
    
  public:
    LagBasis(int n, double a);
    ~LagBasis(){
      delete [] weights;
      delete [] xi;
    }
    
    double phi(int i, double r) const;
    double x(int i)const{ return xi[i];}
    double w(int i)const{ return weights[i];}
};
