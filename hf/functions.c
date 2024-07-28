#include "functions.h"
#include <stdio.h>

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1e-30 //has to be small but not FLP_MIN

void gcf(float *gammcf, float a, float x, float *gln){
    int i;
    float an,b,c,d,del,h;
    *gln=lgammaf(a);
    b=x+1.0-a;
    c=1.0/FPMIN;
    d=1.0/b;
    h=d;

    for (i=1;i<=ITMAX;i++) {
        an = -i*(i-a); b += 2.0;
        d=an*d+b;
        if (fabs(d) < FPMIN)
            d=FPMIN;
        c=b+an/c; if (fabs(c) < FPMIN)
            c=FPMIN;
        d=1.0/d; del=d*c; h *= del;
        if (fabs(del-1.0) < EPS)
            break;
    }
    if (i > ITMAX)
        printf("a too large, ITMAX too small in gcf");

    *gammcf=exp(-x+a*log(x)-(*gln))*h;
}

void gser(float *gamser, float a, float x, float *gln){
    int n; float sum,del,ap;
    *gln=lgammaf(a);
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
        ++ap;
        del *= x/ap;
        sum += del; 
        if (fabs(del) < fabs(sum)*EPS) {
            *gamser=sum*exp(-x+a*log(x)-(*gln));
            return; 
        } 
    } 
    printf("a too large, ITMAX too small in routine gser");
    return;
}

float  gamma_inc(float a, float x){
    //From Numerical recipes in C 2nd ed http://s3.amazonaws.com/nrbook.com/book_C210.html
    //Its the gamma function referred as Q(a,x) in numerical recipes not P(a,x)
    if (x < 0.0 || a <= 0.0) printf("Invalid arguments in routine gammq");
    float gamser,gammcf,gln;
    if (x < (a+1.0)) { 
        gser(&gamser,a,x,&gln);
        return gamser;
    }
    gcf(&gammcf,a,x,&gln);
    return 1.0-gammcf;
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
    if ((double)x < 1e-7){
        return (1.0/(2.0 * fact + 1.0)) - x * (1.0/(2.0*fact+3.0));
    }
    float ft =fact+0.5;
    return 0.5 * powf(x,-ft) * tgammaf(ft) * gamma_inc(ft,x);
}