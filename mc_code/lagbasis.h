
class LagBasis{
  private:
    const int size;
    double * weights;
    double * xi; //zeros of Pn(2x-1)
    
  public:
    LagBasis(int n);
    ~LagBasis(){
      delete [] weights;
      delete [] xi;
    }
    
    double phi(int i, double r);
    double x(int i)const{ return xi[i];}
    double w(int i)const{ return weights[i];}
};
