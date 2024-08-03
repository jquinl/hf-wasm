#include "functions.h"
#include <stdio.h>
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1e-30 //has to be small but not FLP_MIN

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

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

void matrix_dot(const float * M, const float * M2, const int N, float * Mout){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            float a = 0.0f;
            for (int k = 0; k < N; k++){
                a+= M[k + i*N] * M2[j + k*N];
            }
            Mout[j + i*N] = a;
        }
    }
    
}

int lu_decomp_inplace(float * A, int N, int * P){
    //returns int, 0 failed 1 pair number of row exchanges -1 odd number of row exc.
    //Numerical recipes in C
    float V[N];
    int d = 1;
    for (int i=0 ; i<N; i++) {
        float big=0.0;
        for (int j=0; j<N; j++){
            float temp = fabs(A[j+i*N]);
            if (temp>big)
                big=temp;
        }
        if (big == 0.0)
            return 0;
        V[i]=1.0/big;
    }
    for (int j=0; j<N; j++) {
        for (int i=0; i<j; i++) {
            float sum = 0.0f;
            for (int k = 0; k<i; k++)
                sum += A[k+i*N] * A[j+k*N];
            A[j+i*N] -=sum;
        }
        float big =0.0f;
        int imax = 0;
        for (int i = j; i<N; i++) {
            float sum = 0.0f;
            for (int k = 0; k<j; k++) {
                sum += A[k + i*N] * A[j + k*N];
            }
            A[j + i * N] -= sum;
            float tmp = fabsf(A[i + j * N]) * V[i];
            if(tmp > big){
                big = tmp;
                imax = i;
            }
        }
        if(j != imax){
            //swap rows place row w bigger num as pivot
            for (int k = 0; k < N; k++){
                float tmp = A[k + imax *N];
                A[k + imax * N] = A[k+j*N];
                A[k+j*N] =tmp;
            }
            d*=-1;
            
            float tmp = V[imax];//interchage here too
            V[imax] = V[j];
            V[j] = tmp;
        }

        P[j] = imax;
        if(A[j+j*N] == 0) return 0;
        if(j != N-1){
            float tmp = 1.0f/A[j +j * N];
            for(int i=j+1; i< N ; i++){
                A[j+i*N] *= tmp;
            }
        }
    }
    return d;
}

void lu_bcksubs(float * LU,int * IDX,int N ,float * B){

    for (int i = 0; i < N; i++){
        int ip = IDX[i];
        float sum = B[ip];
        B[ip] = B[i];
        for (int  j = 0; j < i; j++){
            sum -= LU[j+i*N] * B[j];
        }
        B[i] = sum;
    }
    for (int i = N-1; i > -1; i--){
        float sum = B[i];
        for (int j = i+1; j < N; j++){
            sum-=LU[j+i*N] * B[j];
        }
        B[i] = sum/LU[i+i*N];
    }
    
}
void matrix_invert(const float * M, const int N, float * Minv){
    float LU[N*N];
    int  IDX[N];
    //memcpy(LU,M,N*N*sizeof(float));
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            LU[j + i *N] = M[j + i*N];
        }
    }

    if(lu_decomp_inplace(LU,N,IDX) == 0){
        printf("Error in LU Decomposition");
    }
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            Minv[j+i*N] = 0.0f;
        }
    }
    float * Ms = Minv;
    for (int i = 0 ; i<N; i++) {
        Minv[i + i * N] = 1.0f;
        lu_bcksubs(LU, IDX, N, Ms);
        Ms+=N;
    }
}

void householder_reduction(float * A, int N, float * D,float * E){
    for (int i = N-1; i > 0; i--){
        int l = i-1;
        float h =0.0f, scale = 0.0f;

        if(l > 0){
            for (int k = 0; k < i; k++){
                scale+=fabs(A[k+i*N]);
            }
            if(scale == 0.0){
                E[i]= A[l+i*N];
            }else{
                for (int k = 0; k < i; k++){
                    A[k+i*N] /=scale;
                    h+=A[k+i*N]*A[k+i*N];
                }
                float f = A[l+i*N];
                float g = (f>= 0.0? -sqrtf(h) : sqrtf(h));
                E[i] = scale * g;
                h -= f*g;
                A[l+i*N] = f-g;
                f =0.0f;
                for (int j = 0; j < i; j++){
                    A[i+j*N] = A[j+i*N]/h;
                    g = 0.0f;
                    for (int k = 0; k < j+1; k++){
                        g+= A[k+j*N] * A[k+i*N];
                    }
                    for (int k = j+1; k < i; k++){
                        g+= A[j+k*N] * A[k+i*N];
                    }
                    E[j] = g/h;
                    f+=E[j]* A[j+i*N];
                }
                float hh = f/(h+h);
                for (int j = 0; j < i; j++){
                    f = A[j+i*N];
                    E[j]=g=E[j]-hh*f;
                    for (int k = 0; k < j+1; k++){
                        A[k+j*N] -= f*E[k]+g*A[k+i*N]; 
                    }
                }
            }
        
        }else{
            E[i] = A[l+i*N];
        }
        D[i] = h;
    }
    D[0]=0.0f;
    E[0]=0.0f;
    for (int i = 0; i < N; i++){
        //int l = i-1 ¯\_(ツ)_/¯ ;
        if(D[i]!=0.0){
            for (int j = 0; j < i; j++){
                float g = 0.0f;
                for (int k = 0; k < i; k++){
                    g+=A[k+i*N] * A[j+k*N];
                }
                for (int k = 0; k < i; k++){
                    A[j+k*N] -= g* A[i+k*N];
                }
            }
        }
        D[i] = A[i+i*N];
        A[i+i*N] = 1.0f;
        for (int j = 0; j < i; j++){
            A[i+j*N] = A[j+i*N] = 0.0f;
        }
    }
}

float pythag(float a, float b)
{
	float absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

int tridiagonal_qli(float * D, float * E, int N, float * Z){
    //returns 1 if fails
    int m,l;
	for (int i=1; i < N; i++){
        E[i-1]=E[i];
    }
    E[N-1] =0.0f;
    for (l = 0; l < N; l++){
        int iter = 0;
        do{
            for (m = l; m < N-1; m++){
               //This seems kindof unstable 
               //seems to more or less scale the vecs equally so ok, some small values get scaled a lot but not the full eigenvec :(
               float dd=fabs(D[m])+fabs(D[m+1]);
               if ((float)(fabs(E[m])+dd) == dd){
                   break; //omg
               }
               //if(fabsf(E[m]/dd) < 5e-11) break;
            }
            if(m != l){
                if(iter++ == 30) return 1;
                float g = (D[l+1]-D[l])/(2.0f*E[l]);
                float r = pythag(g,1.0f);
                g = D[m]-D[l]+E[l] / (g+SIGN(r,g));
                float s=1.0f,c=1.0f;
                float p=0.0f;
                int i = m-1;
                for (; i > l-1; i--){
                    float f = s * E[i];
                    float b = c * E[i];
                    r = pythag(f,g);
                    E[i+1] = r;
                    if(r == 0.0){
                        D[i+1] -=p;
                        E[m] = 0.0f;
                        break;
                    }
                    s = f/r;
                    c = g/r;
                    g = D[i+1]-p;
                    r = (D[i]-g) * s + 2.0f * c * b;
                    p= s*r;
                    D[i+1] = g+p;
                    g = c*r-b;
                    for (int k = 0; k < N; k++){
                        f = Z[i+1+ k*N];
                        Z[i+1 + k*N] = s*Z[i + k*N] + c * f;
                        Z[i + k*N] = c*Z[i + k*N] - s*f;
                    }
                }
                if(r == 0.0 && i>l-1) continue;
                D[l] -=p;
                E[l] = g;
                E[m] = 0.0f;
            }
        }while (m != l);
    }
    return 0;
}

void matrix_eigval(const float * M, const int N, float * Mevec,float * Meval){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            Mevec[j+i*N] = M[j+i*N];
        }
    }
    float E[N];
    householder_reduction(Mevec,N,Meval,E);
    if(tridiagonal_qli(Meval,E,N,Mevec)){
        printf("Failed producing eigenvalues\n");
    }
    
}

void matrix_diag(const float * M, const int N, float * Mdiag){
    float Meval[N];
    matrix_eigval(M,N,Mdiag,Meval);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if( i == j){
                Mdiag[j+i*N] = Meval[i];
            }
            else{
                Mdiag[j+i*N] = 0.0f;
            }
        }
    }
}

void matrix_sqrt_inplace(float * M,int N){
    float Z[N*N];
    float ZT[N*N];
    float D[N*N];
    float DZT[N*N];
    //M = Z⋅D⋅ZT
    matrix_eigval(M,N,Z,D);
    //copy values in top row to diagonal
    
    for (int i = 0; i < N; i++){
        D[i+i*N] = D[i];
    }
    //Zero the rest
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++)
            if(i!=j) D[j+i*N] = 0.0f;
    }
    //sqrtD
    for (int i = 0; i < N; i++){
        D[i+i*N] = sqrtf(D[i+i*N]);
    }
    //zT
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            ZT[j+i*N] = Z[i+j*N];
        }
    }
    matrix_dot(D,ZT,N,DZT);
    matrix_dot(Z,DZT,N,M);
}