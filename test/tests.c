#include "../hf/primitives.h"
#include "../hf/integrals_naive.h"
#include "../hf/functions.h"
#include <stdio.h>
#include <string.h>
//#include "test1_values.h"
//#include "test2_values.h"
//#include "test3_values.h"
#include "test4_values.h"

#define MAX_NUM_CGTO 20
#define MAX_NUM_PRIM MAX_NUM_CGTO * 3
#define EPSILON 1e-5

prim prims[MAX_NUM_PRIM];
cgto cgtos[MAX_NUM_CGTO];
const int cgto_n = cgto_number;

float S[cgto_n * cgto_n];
float Sinv[cgto_n * cgto_n];
float T[cgto_n * cgto_n];
float V[cgto_n * cgto_n];
float EE[cgto_n * cgto_n * cgto_n * cgto_n];

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
void print_matrix(const float * M, int N){
    for (int i = 0; i < N ; i++){
        for (int j = 0; j < N ; j++){
            printf("%f ",  M[i * N + j]);
        }
        printf("\n");
    }
}
int main(){
    int at_num = test_at_num;
    printf("Atom num: %i nCGTO :%i\n",at_num,cgto_n);
    atoms ats[at_num];
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
    lu_decomp_inplace(LU,cgto_n,IDX);
    printf("-----------LU Decomposition------------------\n");
    if(!compare_matrices(LU,LU_test,cgto_n)){
        printf("OK\n");
    }
   
    printf("-----------InvertMatrix------------------\n");
    matrix_invert(S,cgto_n,Sinv);
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
   
    return 0;
};