#include "../primitives.h"
#include "../integrals_naive.h"
#include <stdio.h>
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
float T[cgto_n * cgto_n];
float V[cgto_n * cgto_n];
float EE[cgto_n * cgto_n * cgto_n * cgto_n];

int main(){
    int at_num = test_at_num;
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
    int fail = 0;
    printf("-----------S------------------\n");
    for (int i = 0; i < cgto_n ; i++){
        for (int j = 0; j < cgto_n ; j++){
          //printf("%f ",  S[i * cgto_n + j]);
           if(fabs(S[j + cgto_n * i] - S_test[j + cgto_n* i])> EPSILON){
                printf("WARNING %f is not %f at %i %i  (%i) \n",S[j + cgto_n* i] , S_test[j + cgto_n* i],
                            i,j,j + cgto_n* i );
                fail = 1;
           }
        }
        //printf("\n");
    }
    if(!fail){
        printf("OK\n");
    }
    fail = 0;
    gto_kinetic(cgtos, cgto_n,T);
    printf("-----------T------------------\n");
    for (int i = 0; i < cgto_n ; i++){
        for (int j = 0; j < cgto_n ; j++){
          // printf("%f ",  T[i * cgto_n + j]);
           if(fabs(T[j + cgto_n* i] - T_test[j + cgto_n* i])> EPSILON){
                printf("WARNING %f is not %f at %i %i  (%i) \n",T[j + cgto_n* i] , T_test[j + cgto_n* i],
                            i,j,j + cgto_n* i );
                fail = 1;
           }
        }
        //printf("\n");
    }
    if(!fail){
        printf("OK\n");
    }
    fail = 0;
    gto_electron_nuclear(cgtos, cgto_n,ats,at_num,V);
    printf("-----------V------------------\n");
    for (int i = 0; i < cgto_n ; i++){
        for (int j = 0; j < cgto_n ; j++){
         //  printf("%f ",  V[i * cgto_n + j]);
           if(fabs(V[j + cgto_n* i] - V_test[j + cgto_n* i])> EPSILON){
                printf("WARNING %f is not %f at %i %i  (%i) \n",V[j + cgto_n* i] , V_test[j + cgto_n* i],
                            i,j,j + cgto_n* i );
                fail = 1;
           }
        }
       // printf("\n");
    }
    if(!fail){
        printf("OK\n");
    }
    fail = 0;
    gto_electron_electron(cgtos,  cgto_n,ats,at_num,EE);
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
    return 0;
};