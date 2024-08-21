    
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../hf/primitives.h"
#include "../hf/integrals_naive.h"
#include "../hf/functions.h"
#include "../hf/rhf.h"
#include "../hf/utils.h"

#define cgto_number 2
#define cgto_n cgto_number
#define test_at_num 2

#define MAX_NUM_CGTO 20
#define MAX_NUM_PRIM MAX_NUM_CGTO * 3
#define EPSILON 1e-5

int main(){
    prim prims[MAX_NUM_PRIM];
    cgto cgtos[MAX_NUM_CGTO];

    float S[cgto_n * cgto_n];
    float Sinv[cgto_n * cgto_n];
    float T[cgto_n * cgto_n];
    float V[cgto_n * cgto_n];
    float EE[cgto_n * cgto_n * cgto_n * cgto_n];

    for (int d = 0.0; d < 1000; d++){

        atoms ats[2] = {
            {1,{0.0f,0.0f,0.0f}},
            {1,{(float)d/50.0f ,0.0f,0.0f}}
        };

        build_cgto(test_at_num,ats,cgtos,prims);
        gto_overlap(cgtos,cgto_n,S);
        gto_kinetic(cgtos,cgto_n,T);
        gto_electron_nuclear(cgtos,cgto_n,ats,test_at_num,V);
        gto_electron_electron(cgtos,cgto_n,EE);

        
        //printf("-----------RHF------------------\n");
        
        int n_elec = 0;
        for (int i = 0; i < test_at_num; i++){
            n_elec += ats[i].type;
        }

        float P[cgto_n * cgto_n];
        float Pt[cgto_n * cgto_n];
        float TMP1[100*100];
        for (int i = 0; i < 100; i++)
            for (int j = 0; j < 100; j++)
            TMP1[i+j*100] = 0.0f;
        float rep = nuclear_repulsion(ats,test_at_num);
        restricted_hartree_fock_init(S,T,V,EE,cgto_n,n_elec,P,Sinv,TMP1);
        float de = 1000.0f;
        float olde = 1000.0f+rep;
        while (fabsf(de)>1e-5){
            float ne = rep + restricted_hartree_fock_step(Sinv,T,V,EE,cgto_n,n_elec,P,TMP1);
            de = ne -olde;
            olde = ne;
        }
        float ee = restricted_hartree_fock_step(Sinv,T,V,EE,cgto_n,n_elec,P,TMP1);
       
        printf("[%f,%f,%f,%f],\n",(float)d/50.0f,ee,rep,rep + ee);
    }
    return 0;

}