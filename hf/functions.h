#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <float.h>
#include <math.h>

float gamma_inc(float a, float x);
int factorial(int n);
int double_factorial(int n);
int binomial(int n, int k);
float boys(int fact,float x);

    
void matrix_dot(const float * M, const float * M2, const int N, float * Mout);
void matrix_dotd(const double * M, const double * M2, const int N, double * Mout);
void matrix_invsqrt(const float * M,const int N,float * Minvsqrt,float * TMP );
void matrix_diag(const float * M, const int N, float * Mdiag,float * TMP );
void matrix_eigval(const float * M, const int N, float * Mevec,float * Meval,float * TMP );

#ifdef MY_TESTS
void matrix_invert(const float * M, const int N, float * Minv,float * TMP ,int * ITEMP);
int lu_decomp_inplace(float * A, int N, int * P,float * V);
void householder_reduction(const float * A,int N,float * D,float * E,float * Z);
void householder_reductiond(const double * A,int N,double * D,double * E,double * Z,double * TMP);
int tridiagonal_ql2(float * D, float* E, int N, float * Z);
#endif
#endif
