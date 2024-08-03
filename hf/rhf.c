#include "rhf.h"

void compute_G(float * EE, float * P, int n_basis, float * G){
    for (int i = 0; i < n_basis; i++){
        for (int j = 0; j < n_basis; j++){
            float g = 0.0f;
            for (int k = 0; k < n_basis; k++){
                for (int l = 0; l < n_basis; l++){
                    float J = EE[l + n_basis * (k + n_basis * (j + n_basis * i))];
                    float K = EE[j + n_basis * (k + n_basis * (l + n_basis * i))];
                    g+= (J -0.5f * K) * P[l + n_basis * k];
                }
            }
            G[j + n_basis * i] = g;
        }
    }
}
void compute_F(float * T, float * V, float * G, int n_basis, float * F){
    for (int i = 0; i < n_basis; i++){
        for (int j = 0; j < n_basis; j++){
            F[j + n_basis * i] = T[j + n_basis * i] + V[j + n_basis * i] + G[j + n_basis * i];
        }
    }
}
void compute_P(float * C, int n_basis, float * P){
    //2 electrons x basis
    for (int i = 0; i < n_basis; i++){
        for (int j = 0; j < n_basis; j++){
            float p = 0.0f;
            for (int a = 0; a < n_basis; a++){
                p += 2.0f * C[a + n_basis * i] * C[j + n_basis * a];
            }
            P[j + n_basis * i] = p;
        }
    }
    
}
float hf_expectation(float* P, float * T, float * V, float * F,int n_basis){
    float e = 0.0f;
    for (int i = 0; i < n_basis; i++){
        for (int j = 0; j < n_basis; j++){
            e+=P[(j + i * n_basis)] * (T[(j + i * n_basis)] + V[(j + i * n_basis)] +F[(j + i * n_basis)]);
        }
    }
    return e * 0.5f;
}
float restricted_hartree_fock_step(float * S,float * T, float * V,float * EE,float * P, int n_basis){
    float G[n_basis * n_basis];
    float F[n_basis * n_basis];
    float Sinvsqr[n_basis * n_basis];

    float FS[n_basis * n_basis];
    float SFS[n_basis * n_basis];
    float Evecs[n_basis * n_basis];
    float Evals[n_basis];
    float C[n_basis * n_basis];
    
    //Two electron contrib
    compute_G(EE,P,n_basis,G);
    //Fock matrix (Hcore + P)
    compute_F(T,V,G,n_basis,F);
    //S^{-1}
    matrix_invert(S,n_basis,Sinvsqr);
    //X=S^{-1/2}}
    matrix_sqrt_inplace(Sinvsqr,n_basis);
    //F'=XFX
    matrix_dot(F,Sinvsqr,n_basis,FS);
    matrix_dot(S,FS,n_basis,SFS);
    //Diag F'
    matrix_eigval(SFS,n_basis,Evecs,Evals);
    //Molec ORbs C=XC'
    matrix_dot(Sinvsqr,Evecs,n_basis,C);
    //Update P P_muv = 2\sum^{N/2}{a}C_{\mu a}C^*_{\nu a}
    compute_P(C,n_basis,P);

    //return SCF energy
    return hf_expectation(P,T,V,F,n_basis);
}