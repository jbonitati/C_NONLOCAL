
double coulomb_potential(int z1,int z2,double R,int n);
double central_potential(int l,double R,double u);

/*Perey-Buck and Woods-Saxon Potential*/

double Vd_potential(double Vd, double R,int n,int z,double rvd,double avd);
double  Wd_potential(double Wd,double R,int n,int z,double rwd,double awd);

double V_sp_potential(double Vso,double Rso,double aso,
	double r,double j,int l,int n,int z);
double W_sp_potential(double Wso,double Rwso,double awso,
	double r,double j,int l, int n, int z);

double  V_potential(double Vv,double R,int n,int z,double rv,double av);
double Wv_potential(double Wv,double R,int n,int z,double rwv,double awv);  

std::complex<double> local_potential(double R,int l,double j,int n,
	int z1,int z2,double Vv,double rv,double av,double Wv,
	double rwv,double awv,double Vd,double rvd,double avd,double Wd,
	double rwd,double awd,double Vso,double Rso,double aso,double Wso,
	double Rwso,double awso);
  
std::complex<double> non_local(double r_minus,double r_plus,int l,
	double jj,int n,int z1,int z2, double beta,double Vv1,double rv1,
	double av1,double Wv1,double rwv1,double awv1,double Vd1,double rvd1,
	double avd1,double Wd1,double rwd1,double awd1,double Vso1,
	double Rso1,double aso1,double Wso1,double Rwso1,double awso1);
