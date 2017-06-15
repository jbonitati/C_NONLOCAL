
double columb_potential(int z1,int z2,double R,int n);
double central_potential(int l,double R,double u);

/*Perey-Buck Potential*/

double Vd_potential(double Vd, double R,int n,int z,double rvd,double avd);
double  Wd_potential(double Wd,double R,int n,int z,double rwd,double awd);

double V_sp_potential(double Vso,double Rso,double aso,
	double r,double j,int l,int n,int z);
double W_sp_potential(double Wso,double Rwso,double awso,
	double r,double j,int l, int n, int z);

double  V_potential(double Vv,double R,int n,int z,double rv,double av);
double Wv_potential(double Wv,double R,int n,int z,double rwv,double awv);  
