#include <cstdio>
#include <iostream>
#include <fenv.h>   // enable floating point trap
#include <fstream>
#define MAXN 1000
#define MAXX0 1000
#include <complex>
#include <cmath>
#include <vector>
#include <math.h>       /* pow */
//#include <array>
//#include <vectors.H>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
//#include <Eigen/Dense>
//#include <arrays2.H>
//#include <vectors.H>
#include <iostream>
#include <fstream>
#include <assert.h>  
#include <math.h>
#include "matrix.h"
extern "C" {
    // LU decomoposition of a general matrix
    void CSPTRI(int* M, int *N,  std::complex < double >* A, int* lda, int* IPIV, int* INFO);
    //CSPTRI( UPLO, N, AP, IPIV,	WORK, INFO )
    // generate inverse of a matrix given its LU decomposition
    void ZSPTRI(int* N,  std::complex < double >* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

void inverse( std::complex < double >* A, int N)
{
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    CSPTRI(&N,&N,A,&N,IPIV,&INFO);
    ZSPTRI(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    delete IPIV;
    delete WORK;
}

int main()
{
// const   complex<double> i(0.0,1.0);  
//     std::complex < double >** A [2*2] = {
//        1,2,
//        3,4
//    };

 //   inverse(A, 2);

  //  printf("%f %f\n", A[0], A[1]);
 //   printf("%f %f\n", A[2], A[3]);
   return 0;

}
