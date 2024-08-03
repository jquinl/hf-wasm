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
void matrix_invert(const float * M, const int N, float * Minv);
void matrix_sqrt_inplace(float * M,const int N);
void matrix_diag(const float * M, const int N, float * Mdiag);
void matrix_eigval(const float * M, const int N, float * Mevec,float * Meval);

#ifdef MY_TESTS
int lu_decomp_inplace(float * A, int N, int * P);
void householder_reduction(float * A, int N, float * D,float * E);
int tridiagonal_qli(float * D, float * E, int N, float * Z);
#endif

#endif