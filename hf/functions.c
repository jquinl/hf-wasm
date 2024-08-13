#include "functions.h"
#include <stdio.h>
#include <float.h>
#include "utils.h"
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
    //From Numerical recipes in C 2nd ed 
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
void matrix_dotd(const double * M, const double * M2, const int N, double * Mout){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            Mout[j + i*N] = 0.0;
            for (int k = 0; k < N; k++){
                Mout[j + i*N]+= M[k + i*N] * M2[j + k*N];
            }
            //Mout[j + i*N] = a;
        }
    }
}
#ifdef MY_TESTS
int lu_decomp_inplace(float * A, int N, int * P,float * V){
    //V[N]
    //returns int, 0 failed 1 pair number of row exchanges -1 odd number of row exc.
    //Numerical recipes in C
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
void matrix_invert(const float * M, const int N, float * Minv,float * TMP,int * ITEMP){
    //TEMP total [N*N+N]
    //ITEMP total [N]
    //LU[N*N];
    //IDX[N];
    float * LU = TMP;
    float * TMP2= TMP + N*N;
    int *IDX  = ITEMP;

    //memcpy(LU,M,N*N*sizeof(float));
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            LU[j + i *N] = M[j + i*N];
        }
    }

    if(lu_decomp_inplace(LU,N,IDX,TMP2) == 0){
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
#endif

void householder_reduction(const float *A, int N, float *D, float *E,float *Z){
    //Description:
    //J. H. Wilkinson. (1962-1963). Householder's method for symmetric matrices. , 4(1), 354–361. doi:10.1007/bf01386332 
    //Implementation:
    //Numerical recipes in C 2nd ed
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
           Z[i+j*N]= A[i+j*N];
        }
    }

	for (int i = N - 1; i > 0; i--) {
		int l = i - 1;
        float h = 0.0f;
        float scale = 0.0f;
		if (l > 0) {
			for (int k = 0; k < i; k++) {
				scale += fabsf(Z[k+i*N]); 
			}
			if (scale <FLT_EPSILON) {
				E[i] = Z[l+i*N]; 
			} else {
				for (int k = 0; k < i; k++) {
					Z[k+i*N] /= scale; 
					h += Z[k+i*N]*Z[k+i*N];
				}
				float f =Z[l+i*N]; 
				float g = (f >= 0.0 ? -sqrtf(h) : sqrtf(h));
				E[i] = scale * g;
				h -= f * g;
				Z[l+i*N] = f - g; 
				f = 0.0f;
				for (int j = 0; j < i; j++) {
                    Z[i+j*N] = Z[j+i*N]/ h;
					g = 0.0f;
					for (int k = 0; k < j + 1; k++)
                        g+=Z[k+j*N] * Z[k+i*N];
					for (int k = j + 1; k < i; k++)
                        g+= Z[j+k*N] * Z[k+i*N]; 
					E[j] = g / h;
					f += E[j]  * Z[j+i*N]; 
				}
				float hh = f / (h + h);
				for (int j = 0; j < i; j++) {
					f = Z[j+i*N];
					E[j] = g = E[j] - hh * f;
					for (int k = 0; k < j + 1; k++) {
                        Z[k+j*N]-=(f * E[k] + g * Z[k+i*N]);
					}
				}
			}
		} else {
			E[i] = Z[l+i*N];
		}
		D[i] = h;
	}

	D[0] = 0.0f;
	E[0] = 0.0f;
	for (int i = 0; i < N; i++) {
		if (fabsf(D[i]) > FLT_EPSILON) {
			for (int j = 0; j < i; j++) {
				float g = 0.0f;
				for (int k = 0; k < i; k++) {
					g += Z[k+i*N] * Z[j+k*N];
				}
				for (int k = 0; k < i; k++) {
					Z[j+k*N] -= g * Z[i+k*N];
				}
			}
		}
		D[i] = Z[i+i*N];
		Z[i+i*N]= 1.0f;
		for (int j = 0; j < i; j++) {
            Z[i+j*N] = Z[j+i*N] = 0.0f;
		}
	}
}

int tridiagonal_ql2(float * D, float* E, int N, float * Z){
    //Description/Impl:
    //Hilary Bowdler; R. S. Martin; C. Reinsch; J. H. Wilkinson. (1968). TheQRandQLalgorithms for symmetric matrices. , 11(4), 293–306. doi:10.1007/bf02166681
    //returns 1 if fails

    const float macheps = FLT_EPSILON;

    float b,f;
	for (int i=1;i<N;i++) E[i-1]=E[i];
    E[N-1]= b = f =0.0;

	for (int l=0;l<N;l++) {
        int j = 0;
        float h = macheps * (fabsf(D[l])+fabsf(E[l]));
        if(b<h) b=h;
        int m = l;
        for (; m < N; m++){
            if(fabsf(E[m])<=b) break;
        }
        if(m != l){
            do{
                if(j++ == 30) return 1;
                //shift
                float p = (D[l+1]-D[l])/ (2.0f*E[l]);
                float r =sqrtf(p*p+1.0f);
                h = D[l]-E[l]/(p<0.0f ? p-r : p+r);
                for (int i = l; i < N; i++){
                    D[i] -=h;
                }
                f+=h;
                //Ql transformation
                p = D[m]; 
                float c = 1.0f;
                float s = 0.0f;
                for (int i = m-1; i > l-1; i--){
                    float g = c*E[i];
                    float h = c*p;
                    if(fabsf(p)>= fabsf(E[i])){
                        c = E[i]/p;
                        r = sqrtf(c*c+1.0f);
                        E[i+1] = s*p*r;
                        s = c/r;
                        c= 1.0f/r;
                    }else{
                        c = p/E[i];
                        r = sqrtf(c*c+1.0f);
                        E[i+1] = s*E[i]*r;
                        s=1.0f/r;
                        c/=r;
                    }
                    p=c*D[i]-s*g;
                    D[i+1]= h+s*(c*g+s*D[i]);
                    for (int k = 0; k < N; k++){
                        h = Z[i+1+k*N];
                        Z[i+1+k*N] = s* Z[i+k*N]+c*h;
                        Z[i+k*N] = c* Z[i+k*N]-s*h;
                    }  
                }
                E[l] = s*p;
                D[l] = c*p;
            }while (fabsf(E[l])>b);
        }
        D[l] +=f; 
    }
    //Sort eigvect/values in ascending order for HF
    for (int i = 0; i < N-1; i++){
        int k = i;
        float p = D[k];
        for (int j = i+1; j < N; j++){
            if(D[j]<p){
                k = j; 
                p = D[j];
            }
        }
        if(k!= i){
            D[k] = D[i];
            D[i] = p;
            for (int j = 0; j < N; j++){
                p = Z[i+j*N];
                Z[i+j*N] = Z[k+j*N];
                Z[k+j*N] = p;
            }
        }
    }
    return 0;
}


void matrix_eigval(const float * M, const int N, float * Mevec,float * Meval,float * TMP ){
 
    float * E = TMP;
    householder_reduction(M,N,Meval,E,Mevec);

    if(tridiagonal_ql2(Meval,E,N,Mevec)){
        printf("Failed producing eigenvalues\n");
    }
}

void matrix_diag(const float * M, const int N, float * Mdiag,float *TMP){
    float * Meval = TMP;
    float * TMPT = TMP+N;
    matrix_eigval(M,N,Mdiag,Meval,TMPT);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < i; j++){
            Mdiag[j+i*N] = 0.0f;
            Mdiag[i+j*N] = 0.0f;
        }
        Mdiag[i+i*N] = Meval[i];
    }
}

void matrix_invsqrt(const float * M,int N,float * Minvsqrt,float * TMP){
    //TMP must be at least [3* N *N]
    float * Z = TMP;
    float * ZT = TMP + N*N;
    float * DZT = TMP + 2*N*N;
    //M = Z⋅D⋅ZT D= Minvsqrt
    matrix_eigval(M,N,Z,Minvsqrt,ZT);//TMP=ZT Here N * N+ 2*N < 3*N*N (N>1)

    
    //copy values in top row to diagonal and 1/sqrtD
    for (int i = 0; i < N; i++){
        Minvsqrt[i+i*N] = 1.f/sqrtf(Minvsqrt[i]);
    }
   
    //Zero the rest
    for (int i = 0; i < N; i++){
        for (int j = 0; j < i; j++){
            Minvsqrt[j+i*N] = 0.0f;
            Minvsqrt[i+j*N] = 0.0f;
        }
    }
    //zT
    for (int i = 0; i < N; i++){
        for (int j = i; j < N; j++){
            ZT[j+i*N] = Z[i+j*N];
            ZT[i+j*N] = Z[j+i*N];
        }
    }

    matrix_dot(Minvsqrt,ZT,N,DZT);
    matrix_dot(Z,DZT,N,Minvsqrt);

}