#include <float.h>
#include <math.h>
#include "functions.h"

int gamma_inc(float a, float x, float * gin){
    //Adapted from https://github.com/aromanro/HartreeFock
    float xam, r, sum, ga, t0;
    int k;
    *gin = 0.0;

    if (a < 0.0 || x < 0) return 0;
    else if (x == 0.0) return 1;

    xam = -x + a * log(x);
    
    if (xam > 600 || a > 160.0) return 1;			
    else if (x <= 1.0 + a) {  
        sum = 1.0 / a;
        r = sum;
        for (k = 1; k <= 80; ++k) {
            r *= x / (a + k);
            sum += r;
            if (fabs(r / sum) <= FLT_EPSILON) break; 
        }
        *gin = exp(xam) * sum;
    }
    else {
        t0 = 0.0;
        for (k = 80; k >= 1; --k) 
            t0 = (k - a) / (1.0 + k / (x + t0));

        float gim = exp(xam) / (x + t0);
        ga = tgamma(a);
        *gin = ga - gim;
    }
    return 1;
}

int factorial(int n) {
    if(n < 2)
        return 1;
    return n * factorial(n - 1);
}
int double_factorial(int n) {
    if(n < 2)
        return 1;
    return n * double_factorial(n - 2);
}
int binomial(int n, int k) {
    return factorial(n) / (factorial(k) * factorial(n - k));
}
float boys(int fact,float x){
    if (x > 0.0){
        float igam;
        float ft =fact+0.5;
        gamma_inc(ft,x,&igam);
        return 0.5 * pow(x,-ft) * tgammaf(ft) * igam;
    }
    return (1.0/(2.0 * fact + 1.0)) - x * (1.0/(2.0*fact+3.0));
}