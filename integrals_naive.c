#include <math.h>
#include "hf_settings.h"
//Naive implementation of the molecular integrals following:
//https://rsc.anu.edu.au/~pgill/papers/045Review.pdf
//https://doi.org/10.1021/acs.jchemed.8b00255 (SI)
//The implementations are mostly found in the 4th haoundout of the SI of https://doi.org/10.1021/acs.jchemed.8b00255
//This is a (more or less) 1:1 implementation of the ones found there without any optimizations

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
//Ck term of the overlap
float ck(int la,int lb, float pax,float pbx,int k){
    float coeff;
    for (int i = 0; i < la+1; i++){
        for (int j = 0; j < lb+1; j++){
            if (i+j == k)
                coeff += binomial(la,i) * binomial(lb,j) * pow(pax,la-i) * pow(pbx,lb-j);
        }
    }
    return coeff;    
}
//angular momentum term of the primitive overlap
float sfnx(float inv_a,float p,
            float a_center,float b_center,
            float alpha,float beta,
            int la, int lb){
    float pax = p - a_center;
    float pbx = p - b_center;
    float sx = 0.0f;
    for (int k = 0; k < (int)((la + lb)/2.0f); k++)
    {
        sx +=ck(la,lb,pax,pbx,2 * k) * doble_factorial(2 * k-1) * pow(0.5 * inv_a,k);
    }
    return sx * sqrt(M_PI * inv_a);
}

//Overlap between two primitives
float prim_overlap(prim * prim1, prim * prim2,int * center1,int * center2){
    float l1 = prim1->l,m1 = prim1->m,n1 = prim1->n;
    float l2 = prim2->l,m2 = prim2->m,n2 = prim2->n;
    float alpha = prim1->alpha, beta = prim2->alpha;

    float inv_alpha = 1 /(alpha + beta);
    float alpha_prod = alpha * beta;
    float pfact = pow(M_PI * inv_alpha,1.5f);
    float distx = center1[0] -center2[0];
    float disty = center1[1] -center2[1];
    float distz = center1[2] -center2[2];
    //gaussian product center
    float px = (alpha *  center1[0] + beta * center2[0]) * (inv_alpha);
    float py = (alpha *  center1[1] + beta * center2[1]) * (inv_alpha);
    float pz = (alpha *  center1[2] + beta * center2[2]) * (inv_alpha);
    
    float sq_dist = distx*distx +disty*disty +distz*distz;
    pfact *= exp(-alpha_prod * sq_dist  / inv_alpha);

    pfact *= sfnx(inv_alpha,px,center1[0],center2[0],alpha, beta, l1, l2);
    pfact *= sfnx(inv_alpha,py,center1[1],center2[2],alpha, beta, m1, m2);
    pfact *= sfnx(inv_alpha,pz,center1[2],center2[2],alpha, beta, n1, n2);
    return pfact;
}
//overlap matrix
void gto_overlap(cgto * cgtos, int cgto_num, float * S){
    for (int i = 0; i < cgto_num; i++){
        for (int j = 0; j < cgto_num; j++){
            float s = 0.0f;
            for (int k = 0; k < PRIMITIVES_PER_ORBITAL; k++){
                for (int l = 0; l < PRIMITIVES_PER_ORBITAL; l++){
                    s += prim_overlap(&(cgtos[i].primitives[k]),&(cgtos[j].primitives[l]),
                                    cgtos[i].center,cgtos[j].center) *
                                    cgtos[i].primitives[k].coeff * cgtos[j].primitives[l].coeff *
                                    cgtos[i].primitives[k].normf * cgtos[j].primitives[l].normf;
                }
            }
            S[i,j] = s;
        }
    }
    
}
//Overlap between two primitives
float prim_kinetic(prim * prim1, prim * prim2,int * center1,int * center2){
    float l1 = prim1->l,m1 = prim1->m,n1 = prim1->n;
    float l2 = prim2->l,m2 = prim2->m,n2 = prim2->n;
    float alpha = prim1->alpha, beta = prim2->alpha;

    float inv_alpha = 1 /(alpha + beta);
    float alpha_prod = alpha * beta;

    float px = (alpha *  center1[0] + beta * center2[0]) * (inv_alpha);
    float py = (alpha *  center1[1] + beta * center2[1]) * (inv_alpha);
    float pz = (alpha *  center1[2] + beta * center2[2]) * (inv_alpha);
    
    float pfact = pow(M_PI * inv_alpha,1.5f);

    //Precalculate some repeated integrals
    float sx = sfnx(inv_alpha, px, center1[0], center2[0], alpha, beta, l1, l2) ;
    float sy = sfnx(inv_alpha, py, center1[0], center2[0], alpha, beta, m1, m2) ;
    float sz = sfnx(inv_alpha, pz, center1[0], center2[0], alpha, beta, n1, n2) ;
    //Term 1
    float k  = beta * (2 * (l2+m2+n2) +3) * sx * sy * sz;
   
    //Term 2
    float beta2 = beta * beta;
    k -= (2 * beta2) * sfnx(inv_alpha,px,center1[0],center2[0],alpha , beta,l1,l2 +2) * sy * sz;
    k -= (2 * beta2) * sx * sfnx(inv_alpha, py, center1[0], center2[0], alpha, beta, m1, m2 +2) * sz;
    k -= (2 * beta2) * sx * sy * sfnx(inv_alpha, pz, center1[0], center2[0], alpha, beta, n1, n2 + 2);

    //Term 3
    k -= 0.5 * (l2 * (l2-1)) * sfnx(inv_alpha,px,center1[0],center2[0],alpha , beta,l1,l2 -2) * sy * sz;
    k -= 0.5 * (m2 * (m2-1)) * sx * sfnx(inv_alpha, py, center1[0], center2[0], alpha, beta, m1, m2 -2) * sz;
    k -= 0.5 * (n2 * (n2-1)) * sx * sy * sfnx(inv_alpha, pz, center1[0], center2[0], alpha, beta, n1, n2 - 2);

    return k;
}

void gto_kinetic(cgto * cgtos, int cgto_num,float * K){
     for (int i = 0; i < cgto_num; i++){
        for (int j = 0; j < cgto_num; j++){
            float k = 0.0f;
            for (int k = 0; k < PRIMITIVES_PER_ORBITAL; k++){
                for (int l = 0; l < PRIMITIVES_PER_ORBITAL; l++){
                    k += prim_kinetic(&(cgtos[i].primitives[k]),&(cgtos[j].primitives[l]),
                                    cgtos[i].center,cgtos[j].center) *
                                    cgtos[i].primitives[k].coeff * cgtos[j].primitives[l].coeff *
                                    cgtos[i].primitives[k].normf * cgtos[j].primitives[l].normf;
                }
            }
            K[i,j] = k;//Leave here
        }
    }
}