#include "rhf.h"
#include <stdio.h>


void compute_F(const float * T, const float * V, const float * EE,const float * P,const int n_basis, float * F){
    for (int i = 0; i < n_basis; i++){
        for (int j = 0; j < n_basis; j++){
            float f = 0.0f;
            for (int k = 0; k < n_basis; k++){
                for (int l = 0; l < n_basis; l++){
                    float J = EE[l + n_basis * (k + n_basis * (j + n_basis * i))];
                    float K = EE[j + n_basis * (k + n_basis * (l + n_basis * i))];
                    f+= (J -0.5f * K) * P[l + n_basis * k];
                }
            }
            F[j + i* n_basis ] = f+ T[j + n_basis * i] + V[j + n_basis * i];
        }
    }
}
void compute_P(const float * C, const int n_basis, const int n_elec, float * P){
    for (int i = 0; i < n_basis; i++){
        for (int j = 0; j < n_basis; j++){
            float p = 0.0f;
            for (int a = 0; a < (int)(n_elec/2); a++){
                p += C[a + i* n_basis ] * C[a + j*n_basis ];
            }
            P[j + i* n_basis ] =  2.0f * p;
        }
    }
    
}
float hf_expectation(const float* P, const float * T, const float * V, const float * F,const int n_basis){
    float e = 0.0f;
    for (int i = 0; i < n_basis; i++){
        for (int j = 0; j < n_basis; j++){
            e+= P[i + j * n_basis] * (T[j + i * n_basis] + V[j + i * n_basis] + F[j + i * n_basis]);
        }
    }
    return e * 0.5f ;
}

void restricted_hartree_fock_init(const float * S,const float * T, const float * V,const float * EE,const int n_basis,const int n_elec,float * P,float* Sinvsqr,float * TMP){
    //X=S^{-1/2}}
    matrix_invsqrt(S,n_basis,Sinvsqr,TMP);
    //zero out P 
    for (int i = 0; i < n_basis; i++){
        for (int j = 0; j < n_basis; j++){
            P[j+i*n_basis] = 0.0f;
        }
    }
    
}
float nuclear_repulsion(atoms * ats, int at_num){
    float r = 0.0f;
    for (int i = 0; i < at_num; i++){
        for (int j = i+1; j < at_num; j++){
            float dx = ats[i].pos[0] -ats[j].pos[0];
            float dy = ats[i].pos[1] -ats[j].pos[1];
            float dz = ats[i].pos[2] -ats[j].pos[2];
            float d = sqrtf(dx*dx+dy*dy+dz*dz);
            r += (float)(ats[i].type * ats[j].type) / d;
        }
    }
    return r;
}
float restricted_hartree_fock_step( const float * Sinvsqr,const  float * T, const  float * V, const float * EE,const int n_basis,const int n_elec,float * P,float * TMP){
    ///Necessary size for TMP 5*N*N+N + [N * N+ 2*N](Eigval) = 6*N*N+3*N
    int nb2 = n_basis*n_basis;
    float* F = TMP;
    float* Fp = F + nb2;
    float* Cp = Fp +  nb2;
    float* C = Cp +  nb2;
    float* Evals = C + nb2 ;
    float* TMPSUB = Evals + n_basis ;

    for (int i = 0; i < nb2; i++)
    {
        F[i] = 0.0f;
        Fp[i] = 0.0f;
        Cp[i] = 0.0f;
        C[i] = 0.0f;
    }
    for (int i = 0; i < n_basis; i++){
        Evals[i] = 0.0f;
    }
    
    
    //print_matrix(P,n_basis);

    //Fock matrix (Hcore+G + P) all in one
    compute_F(T,V,EE,P,n_basis,F);

    //F'=XTFX use TMPSUB as intermediary
    matrix_dot(F,Sinvsqr,n_basis,TMPSUB);
   
    matrix_dot(Sinvsqr,TMPSUB,n_basis,Fp);
 
    //symmetrize
    for (int i = 0; i < n_basis; i++)
    for (int j = 0; j < i; j++){
        Fp[j+i*n_basis] = 0.5*(Fp[j+i*n_basis]+Fp[i+j*n_basis]);
        Fp[j+i*n_basis] = Fp[j+i*n_basis];
    }
   // printf("Fp--------\n");
    //print_matrix(Fp,n_basis);

    //Diag F'C'=EC' SinvT = C' (reused)
    matrix_eigval(Fp,n_basis,Cp,Evals,TMPSUB);

    //printf("Cp--------\n");
   // print_matrix(Cp,n_basis);
    //Molec ORbs C=XC'
    matrix_dot(Sinvsqr,Cp,n_basis,C);

    //Update P 
    compute_P(C,n_basis,n_elec,P);
    for (int i = 0; i < n_basis; i++)
    for (int j = 0; j < i; j++){
        P[j+i*n_basis] = 0.5*(P[j+i*n_basis]+P[i+j*n_basis]);
        P[j+i*n_basis] = P[j+i*n_basis];
    }
    //return SCF energy
    return hf_expectation(P,T,V,F,n_basis);
}