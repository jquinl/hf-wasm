#include <math.h>
#include "functions.h"
#include "hf_settings.h"
//Naive implementation of the molecular integrals following:
//https://rsc.anu.edu.au/~pgill/papers/045Review.pdf
//https://doi.org/10.1021/acs.jchemed.8b00255 (SI)
//The implementations are mostly found in the 4th haoundout of the SI of https://doi.org/10.1021/acs.jchemed.8b00255
//This is a (more or less) 1:1 implementation of the ones found there without any optimizations

//Ck term of the overlap
float ck(int la,int lb, float pax,float pbx,int k){
    float coeff = 0.0;
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
        sx +=ck(la,lb,pax,pbx,2 * k) * double_factorial(2 * k-1) * pow(0.5 * inv_a,k);
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
            for (int k = 0; k < cgtos[i].n_prim; k++){
                for (int l = 0; l < cgtos[j].n_prim; l++){
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

    float px = (alpha *  center1[0] + beta * center2[0]) * (inv_alpha);
    float py = (alpha *  center1[1] + beta * center2[1]) * (inv_alpha);
    float pz = (alpha *  center1[2] + beta * center2[2]) * (inv_alpha);
    
    float pfact = pow(M_PI * inv_alpha,1.5f);

    //Precalculate some repeated integrals
    float sx = sfnx(inv_alpha, px, center1[0], center2[0], alpha, beta, l1, l2) ;
    float sy = sfnx(inv_alpha, py, center1[0], center2[0], alpha, beta, m1, m2) ;
    float sz = sfnx(inv_alpha, pz, center1[0], center2[0], alpha, beta, n1, n2) ;
    //Term 1
    float t  = beta * (2 * (l2+m2+n2) +3) * sx * sy * sz;
   
    //Term 2
    float beta2 = beta * beta;
    t -= (2 * beta2) * sfnx(inv_alpha,px,center1[0],center2[0],alpha , beta,l1,l2 +2) * sy * sz;
    t -= (2 * beta2) * sx * sfnx(inv_alpha, py, center1[0], center2[0], alpha, beta, m1, m2 +2) * sz;
    t -= (2 * beta2) * sx * sy * sfnx(inv_alpha, pz, center1[0], center2[0], alpha, beta, n1, n2 + 2);

    //Term 3
    t -= 0.5 * (l2 * (l2-1)) * sfnx(inv_alpha,px,center1[0],center2[0],alpha , beta,l1,l2 -2) * sy * sz;
    t -= 0.5 * (m2 * (m2-1)) * sx * sfnx(inv_alpha, py, center1[0], center2[0], alpha, beta, m1, m2 -2) * sz;
    t -= 0.5 * (n2 * (n2-1)) * sx * sy * sfnx(inv_alpha, pz, center1[0], center2[0], alpha, beta, n1, n2 - 2);

    return pfact * t;
}

void gto_kinetic(cgto * cgtos, int cgto_num,float * T){
    for (int i = 0; i < cgto_num; i++){
        for (int j = 0; j < cgto_num; j++){
            float t = 0.0f;
            for (int k = 0; k < cgtos[i].n_prim; k++){
                for (int l = 0; l < cgtos[j].n_prim; l++){
                    t += prim_kinetic(&(cgtos[i].primitives[k]),&(cgtos[j].primitives[l]),
                                    cgtos[i].center,cgtos[j].center) *
                                    cgtos[i].primitives[k].coeff * cgtos[j].primitives[l].coeff *
                                    cgtos[i].primitives[k].normf * cgtos[j].primitives[l].normf;
                }
            }
            T[i,j] = t;//Leave here
        }
    }
}

float vlri(int l, int r, int i, int la,int lb, float  PAx, float  PBx,float PCx, float eps){
    float vlri = pow(-1,l);
    vlri *= ck(la,lb,PAx,PBx,l);
    vlri *= pow(-1,i) * factorial(l) * pow(PCx,l-2*r-2*i) * pow(eps,r+i);
    vlri /=  factorial(r) *factorial(i) * factorial(l-2*r-2*i);
    return vlri;
}
float prim_electron_nuclear(prim * prim1, prim * prim2,int * center1,int * center2,int Z, int * Zr){
    float l1 = prim1->l,m1 = prim1->m,n1 = prim1->n;
    float l2 = prim2->l,m2 = prim2->m,n2 = prim2->n;
    float alpha = prim1->alpha, beta = prim2->alpha;

    float inv_alpha = 1 /(alpha + beta);
    float alpha_prod = alpha * beta;
    //gaussian product center
    float px = (alpha *  center1[0] + beta * center2[0]) * (inv_alpha);
    float py = (alpha *  center1[1] + beta * center2[1]) * (inv_alpha);
    float pz = (alpha *  center1[2] + beta * center2[2]) * (inv_alpha);
    //squared distance AB
    float distx = center1[0] -center2[0];
    float disty = center1[1] -center2[1];
    float distz = center1[2] -center2[2];
    float sq_dist = distx*distx +disty*disty +distz*distz;
    //squared distance PC
    float pcx = px - Zr[0];
    float pcy = py - Zr[1];
    float pcz = pz - Zr[2];
    float sqPC = pcx*pcx +pcy*pcy +pcz*pcz;
    //PA,PB
    float pax = px-center1[0];
    float pay = py-center1[1];
    float paz = pz-center1[2];
    float pbx = px-center2[0];
    float pby = py-center2[1];
    float pbz = pz-center2[2];
    //Term for the boys fn
    float fterm = (alpha + beta) * sqPC;

    float v = -Z * (2* M_PI * inv_alpha  ) * exp(-alpha_prod * sq_dist * inv_alpha) ;
    float eps = inv_alpha * 0.25;
    float v_accum = 0.0f;
    for (int l = 0; l < l1  +l2; l++){
        for (int r = 0; r < l/2; r++){
            for (int i = 0; i < (l-2* r)/2; i++){
                float vx = vlri(l,r,i,l1,l2,pax,pbx,pcx,eps);
                for (int m = 0; m < m1  +m2; m++){
                    for (int s = 0; s < m/2; s++){
                        for (int j = 0; j < (m-2* s)/2; j++){
                            float vy = vlri(m,s,j,m1,m2,pay,pby,pcy,eps);
                            for (int n = 0; n < n1  +n2; n++){
                                for (int t = 0; t < n/2; t++){
                                    for (int k = 0; k < (n-2* t)/2; k++){
                                        float vz = vlri(n,t,k,n1,n2,paz,pbz,pcz,eps);
                                        int fsub = l+m+n -2 * (r+s+t)-(i+j+k);
                                        float F = boys(fsub,fterm);
                                        v_accum += vx * vy * vz * F;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    v *= v_accum;
    return v;
}

void gto_electron_nuclear(cgto * cgtos, int cgto_num, atoms * atom_array, int atom_num,float * V){
        
    for (int i = 0; i < cgto_num; i++){
        for (int j = 0; j < cgto_num; j++){
            for (int  a = 0; a < atom_num; a++)
            {
                float v = 0.0f;
                for (int k = 0; k < cgtos[i].n_prim; k++){
                    for (int l = 0; l < cgtos[j].n_prim; l++){
                        v += prim_electron_nuclear(&(cgtos[i].primitives[k]),&(cgtos[j].primitives[l]),
                                        cgtos[i].center,cgtos[j].center, atom_array[a].type, atom_array[a].pos) *
                                        cgtos[i].primitives[k].coeff * cgtos[j].primitives[l].coeff *
                                        cgtos[i].primitives[k].normf * cgtos[j].primitives[l].normf;
                    }
                }
                V[i,j] += v;//Leave here
            }
        }
    }
}