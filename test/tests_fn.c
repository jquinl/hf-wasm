#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../hf/primitives.h"
#include "../hf/integrals_naive.h"
#include "../hf/functions.h"
#include "../hf/rhf.h"
#include "../hf/utils.h"

#include "test_values/test1_values.h"
//#include "test_values/test2_values.h"
//#include "test_values/test3_values.h"
//#include "test_values/test4_values.h"

#define MAX_NUM_CGTO 20
#define MAX_NUM_PRIM MAX_NUM_CGTO * 3
#define EPSILON 1e-5

prim prims[MAX_NUM_PRIM];
cgto cgtos[MAX_NUM_CGTO];

float S[cgto_n * cgto_n];
float Sinv[cgto_n * cgto_n];
float T[cgto_n * cgto_n];
float V[cgto_n * cgto_n];
float EE[cgto_n * cgto_n * cgto_n * cgto_n];

float TMP[7 * cgto_n * cgto_n+3*cgto_n];
int ITMP[cgto_n ];


int compare_matrices(const float * M1,const float * M2, int N){
    int fail = 0;
    for (int i = 0; i < N ; i++){
        for (int j = 0; j < N ; j++){
            if(fabs(M1[j + N * i] - M2[j + N* i])> EPSILON){
                printf("WARNING %f is not %f at %i %i  (%i) \n",M1[j + N* i] , M2[j + N* i],
                            i,j,j + N* i );
                fail = 1;
            }
        }
    }
    return fail;
}
int compare_matrices_abs(const float * M1,const float * M2, int N){
    int fail = 0;
    for (int i = 0; i < N ; i++){
        for (int j = 0; j < N ; j++){
            if(fabsf(fabsf(M1[j + N * i]) - fabsf(M2[j + N* i]))> EPSILON){
                printf("WARNING %f is not %f at %i %i  (%i) \n",M1[j + N* i] , M2[j + N* i],
                            i,j,j + N* i );
                fail = 1;
            }
        }
    }
    return fail;
}
int compare_matricesd(const double * M1,const double * M2, int N){
    int fail = 0;
    for (int i = 0; i < N ; i++){
        for (int j = 0; j < N ; j++){
            if(fabs(M1[j + N * i] - M2[j + N* i])> EPSILON){
                printf("WARNING %f is not %f at %i %i  (%i) \n",M1[j + N* i] , M2[j + N* i],
                            i,j,j + N* i );
                fail = 1;
            }
        }
    }
    return fail;
}
int compare_vectors(const float * V1,const float * V2, int N){
    int fail = 0;
    for (int i = 0; i < N ; i++){
        if(fabs(V1[i] - V2[i])> EPSILON){
            printf("WARNING %f is not %f at %i \n",V1[i] , V2[i],i );
            fail = 1;
        }
    }
    return fail;
}
int compare_vectorsd(const double * V1,const double * V2, int N){
    int fail = 0;
    for (int i = 0; i < N ; i++){
        if(fabs(V1[i] - V2[i])> EPSILON){
            printf("WARNING %f is not %f at %i \n",V1[i] , V2[i],i );
            fail = 1;
        }
    }
    return fail;
}

void make_simmetric(float * M, int N){
    for (int i = 0; i < N; i++){
        for (int j = i; j < N; j++){
            M[j+i*N] = 0.5f * ( M[j+i*N] + M[i+j*N] );
        }
    }
    
}
int main(){
    int at_num = test_at_num;
    printf("Atom num: %i nCGTO :%i\n",at_num,cgto_n);
    atoms ats[test_at_num];
    for (int i = 0; i < at_num; i++){
        ats[i] = (atoms){test_types[i],{test_pos[i][0],test_pos[i][1],test_pos[i][2]}};
    }


    build_cgto(at_num,ats,cgtos,prims);
   //printf("-------------\n");
   //for (int i = 0; i < 6 * ORBITALS_PER_ATOM; i++)
   //{
   //    printf("%f \n", prims[i].alpha);
   //    
   //}
    gto_overlap(cgtos, cgto_n,S);
    printf("-----------S------------------\n");
    if(!compare_matrices(S,S_test,cgto_n)){
        printf("OK\n");
    }
    //print_matrix(S,cgto_n);


    gto_kinetic(cgtos, cgto_n,T);
    printf("-----------T------------------\n");
    if(!compare_matrices(T,T_test,cgto_n)){
        printf("OK\n");
    }
    gto_electron_nuclear(cgtos, cgto_n,ats,at_num,V);
    printf("-----------V------------------\n");
    if(!compare_matrices(V,V_test,cgto_n)){
        printf("OK\n");
    }
    int fail = 0;
    gto_electron_electron(cgtos,  cgto_n,EE);
    printf("-----------EE------------------\n");
    for (int i = 0; i < cgto_n ; i++){
        for (int j = 0; j < cgto_n ; j++){
           // printf("-----------------\n");
            for (int k = 0; k < cgto_n ; k++){
                for (int l = 0; l < cgto_n ; l++){
                   //printf("%f ",  EE[l + cgto_n * (k + cgto_n * (j + cgto_n* i))]);
                   if(fabs(EE[l + cgto_n * (k + cgto_n * (j + cgto_n* i))] - EE_test[l + cgto_n * (k + cgto_n * (j + cgto_n* i))])> EPSILON){
                        printf("WARNING %f is not %f at %i %i %i %i (%i) ",EE[l + cgto_n * (k + cgto_n * (j + cgto_n* i))] , EE_test[l + cgto_n * (k + cgto_n * (j + cgto_n* i))],
                            i,j,k,l,l + cgto_n * (k + cgto_n * (j + cgto_n* i)) );
                        fail = 1;
                   }
                }
            //printf("\n");
            }
        }
       // printf("\n");
    }
    if(!fail){
        printf("OK\n");
    }

    float LU[cgto_n * cgto_n];
    memcpy(LU, S,cgto_n * cgto_n * sizeof(float));
    int IDX[cgto_n ];
    lu_decomp_inplace(LU,cgto_n,IDX,TMP);
    printf("-----------LU Decomposition------------------\n");
    if(!compare_matrices(LU,LU_test,cgto_n)){
        printf("OK\n");
    }

    printf("-----------InvertMatrix------------------\n");
    matrix_invert(S_test,cgto_n,Sinv,TMP,ITMP);
    if(!compare_matrices(Sinv,INV_test,cgto_n)){
        printf("OK\n");
    }
    printf("-----------Matmul------------------\n");
    float ID_test[cgto_n * cgto_n] ;
    for (int i = 0; i < cgto_n ; i++){
        for (int j = 0; j < cgto_n ; j++){
            ID_test[i * cgto_n + j] = 0.0f;
            if (i== j){
                ID_test[i * cgto_n + j] = 1.0f;
            }
        }
    }
    float ID[cgto_n * cgto_n];
    matrix_dot(S,Sinv,cgto_n,ID);
    if(!compare_matrices(ID,ID_test,cgto_n)){
        printf("OK\n");
    }

    printf("-----------Eigval------------------\n");
    float Evec[cgto_n*cgto_n];
    float Evals[cgto_n];
    int failed = 0;
    make_simmetric(S_test,cgto_n);
    
    matrix_eigval(S_test, cgto_n, Evec,Evals,TMP);
    for (int i = 0; i < cgto_n; i++){
        failed = 1;
        if(fabsf(Evals[i] - Eval_test[i])< EPSILON){
            failed =0;
        }
        if(failed){
            break;
        }
    }
    if(!failed){
        printf("OK\n");
    }
    printf("-----------Eigvec------------------\n");
    if(!compare_matrices_abs(Evec,Evec_test,cgto_n)){
        printf("OK\n");
    }

    printf("-----------MSqrt------------------\n");
    float Ssqrt[cgto_n *cgto_n];
    float Sout[cgto_n *cgto_n];
    float Sinv[cgto_n *cgto_n];
    matrix_invert(S_test,cgto_n,Sinv,TMP,ITMP);
    matrix_invsqrt(S_test,cgto_n,Ssqrt,TMP);
    matrix_dot(Ssqrt,Ssqrt,cgto_n,Sout);
    if(!compare_matrices(Sout,Sinv,cgto_n)){
        printf("OK\n");
    }


    return 0;
};