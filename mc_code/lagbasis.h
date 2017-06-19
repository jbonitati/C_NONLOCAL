
class LagBasis{
  private:
    double * weights;
    double * xi; //zeros of Pn(2x-1)
    
  public:
    LagBasis(int n);
    ~LagBasis(){
      delete [] weights;
      delete [] xi;
    }
    
};
