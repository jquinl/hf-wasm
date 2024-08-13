#include "utils.h"

void print_vector(const float * V, int N){
    for (int i = 0; i < N ; i++){
        printf("%f ",  V[i]);
    }
    printf("\n");
}
void print_vectord(const double * V, int N){
    for (int i = 0; i < N ; i++){
        printf("%f ",  V[i]);
    }
    printf("\n");
}
void print_matrix(const float * M, int N){
    for (int i = 0; i < N ; i++){
        for (int j = 0; j < N ; j++){
            printf("%f ",  M[i * N + j]);
        }
        printf("\n");
    }
}
void print_matrixd(const double * M, int N){
    for (int i = 0; i < N ; i++){
        for (int j = 0; j < N ; j++){
            printf("%f ",  M[i * N + j]);
        }
        printf("\n");
    }
}