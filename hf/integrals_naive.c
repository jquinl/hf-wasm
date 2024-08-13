
#include <math.h>
#include "integrals_naive.h"
#include <stdio.h>
#include <stdlib.h>
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
                coeff += binomial(la,i) * binomial(lb,j) * powf(pax,abs(la-i)) * powf(pbx,abs(lb-j)); //abs so it doesnt nan
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
    for (int k = 0; k < (int)((la + lb)/2.0f)+1; k++)
    {
        sx +=ck(la,lb,pax,pbx,2 * k) * double_factorial(2 * k-1) * powf(0.5 * inv_a,k);
    }
    return sx * sqrt(M_PI * inv_a);
}

//Overlap between two primitives
float prim_overlap(prim * prim1, prim * prim2, float * center1, float * center2){
    int l1 = prim1->l,m1 = prim1->m,n1 = prim1->n;
    int l2 = prim2->l,m2 = prim2->m,n2 = prim2->n;
    float alpha = prim1->alpha, beta = prim2->alpha;

    float inv_alpha = 1.0f /(alpha + beta);
    float alpha_prod = alpha * beta;
    float distx = center1[0] - center2[0];
    float disty = center1[1] - center2[1];
    float distz = center1[2] - center2[2];
    //gaussian product center
    float px = (alpha *  center1[0] + beta * center2[0]) * (inv_alpha);
    float py = (alpha *  center1[1] + beta * center2[1]) * (inv_alpha);
    float pz = (alpha *  center1[2] + beta * center2[2]) * (inv_alpha);
    float sq_dist = distx*distx +disty*disty +distz*distz;
    float pfact= exp(-alpha_prod * sq_dist  * inv_alpha);
    pfact *= sfnx(inv_alpha,px,center1[0],center2[0],alpha, beta, l1, l2);
    pfact *= sfnx(inv_alpha,py,center1[1],center2[1],alpha, beta, m1, m2);
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
            S[j + cgto_num * i] = s;
        }
    }
    
}
//Overlap between two primitives
float prim_kinetic(prim * prim1, prim * prim2,float * center1,float * center2){
    int l1 = prim1->l,m1 = prim1->m,n1 = prim1->n;
    int l2 = prim2->l,m2 = prim2->m,n2 = prim2->n;
    float alpha = prim1->alpha, beta = prim2->alpha;

    float inv_alpha = 1.0f /(alpha + beta);
    float alpha_prod = alpha * beta;
    float distx = center1[0] - center2[0];
    float disty = center1[1] - center2[1];
    float distz = center1[2] - center2[2];
    float px = (alpha *  center1[0] + beta * center2[0]) * (inv_alpha);
    float py = (alpha *  center1[1] + beta * center2[1]) * (inv_alpha);
    float pz = (alpha *  center1[2] + beta * center2[2]) * (inv_alpha);
    float sq_dist = distx*distx +disty*disty +distz*distz;
    
    float pfact= exp(-alpha_prod * sq_dist  * inv_alpha);

    //Precalculate some repeated integrals
    float sx = sfnx(inv_alpha, px, center1[0], center2[0], alpha, beta, l1, l2) ;
    float sy = sfnx(inv_alpha, py, center1[1], center2[1], alpha, beta, m1, m2) ;
    float sz = sfnx(inv_alpha, pz, center1[2], center2[2], alpha, beta, n1, n2) ;
    //Term 1
    float t  = beta * (2 * (l2+m2+n2) +3) * sx * sy * sz;
   
    //Term 2
    float beta2 = beta * beta;
    t -= (2 * beta2) * sfnx(inv_alpha,px,center1[0],center2[0],alpha , beta,l1,l2 +2) * sy * sz;
    t -= (2 * beta2) * sx * sfnx(inv_alpha, py, center1[1], center2[1], alpha, beta, m1, m2 +2) * sz;
    t -= (2 * beta2) * sx * sy * sfnx(inv_alpha, pz, center1[2], center2[2], alpha, beta, n1, n2 + 2);

    //Term 3
    t -= 0.5 * (l2 * (l2-1)) * sfnx(inv_alpha,px,center1[0],center2[0],alpha , beta,l1,l2 -2) * sy * sz;
    t -= 0.5 * (m2 * (m2-1)) * sx * sfnx(inv_alpha, py, center1[1], center2[1], alpha, beta, m1, m2 -2) * sz;
    t -= 0.5 * (n2 * (n2-1)) * sx * sy * sfnx(inv_alpha, pz, center1[2], center2[2], alpha, beta, n1, n2 - 2);

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
            T[j + cgto_num * i] = t;//Leave here
        }
    }
}

float vlri(int l, int r, int i, int la,int lb, float  PAx, float  PBx,float PCx, float eps){
    float vlri = powf(-1,l);
    vlri *= ck(la,lb,PAx,PBx,l);
    vlri *= powf(-1,i) * factorial(l) * powf(PCx,l-2*r-2*i) * powf(eps,r+i);
    vlri /=  factorial(r) *factorial(i) * factorial(l-2*r-2*i);
    return vlri;
}
float prim_electron_nuclear(prim * prim1, prim * prim2,float * center1,float * center2,int Z, float * Zr){
    int l1 = prim1->l,m1 = prim1->m,n1 = prim1->n;
    int l2 = prim2->l,m2 = prim2->m,n2 = prim2->n;
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
    float fterm = fabsf((alpha + beta) * sqPC);
    
    
    float v = -Z * (2* M_PI * inv_alpha  ) * exp(-alpha_prod * sq_dist * inv_alpha) ;
    float eps = inv_alpha * 0.25;
    float v_accum = 0.0f;
    for (int l = 0; l < l1+l2 +1; l++){
        for (int r = 0; r < (int)(l/2) +1; r++){
            for (int i = 0; i < (int)((l-2* r)/2) +1; i++){
                float vx = vlri(l,r,i,l1,l2,pax,pbx,pcx,eps);

                for (int m = 0; m < m1+m2 +1; m++){
                    for (int s = 0; s < (int)(m/2) +1; s++){
                        for (int j = 0; j < (int)((m-2* s)/2) +1; j++){
                            float vy = vlri(m,s,j,m1,m2,pay,pby,pcy,eps);
                            for (int n = 0; n < n1+n2 +1; n++){
                                for (int t = 0; t < (int)(n/2) +1; t++){
                                    for (int k = 0; k < (int)((n-2* t)/2) + 1; k++){
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
            V[j + cgto_num * i] = 0.0f;
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
                V[j + cgto_num * i] += v;//Leave here
            }
        }
    }
}

float theta(int l, int la, int lb, float a, float b,int r, float alpha_sum){

    return ck(la,lb,a,b,l) * factorial(l) * powf(alpha_sum,r-l) / (factorial(r)* factorial(l-2*r));
}
float gfn(int lp, int lq,int rp,int rq,int i,
            int la,int lb, int lc,int ld,
            float pa,float pb,float qc, float qd,
            float pq,float inv1,float inv2,float delta){
    float g = powf(-1,i) * powf(-1,lp) * powf(2*delta,2*(rp+rq)) * factorial(lp+lq - 2* rp -2*rq)* powf(delta,i) *
                powf(pq,lp+lq-2* (rp+rq+i));
    g /= powf(4*delta,lp+lq) * factorial(i) * factorial(lp+lq-2*(rp+rq+i));
    g *= theta(lp,la,lb,pa,pb,rp,1/inv1) * theta(lq,lc,ld,qc,qd,rq,1/inv2);
    return g;
}
float  prim_electron_ee(prim * prim1, prim * prim2,prim * prim3, prim * prim4,
                            float * center1, float * center2, float * center3, float * center4){
    int l1 = prim1->l,m1 = prim1->m,n1 = prim1->n;
    int l2 = prim2->l,m2 = prim2->m,n2 = prim2->n;
    int l3 = prim3->l,m3 = prim3->m,n3 = prim3->n;
    int l4 = prim4->l,m4 = prim4->m,n4 = prim4->n;
    float alpha1 = prim1->alpha, alpha2 = prim2->alpha;
    float alpha3 = prim3->alpha, alpha4 = prim4->alpha;
    float inv_alpha1 = 1.0f /(alpha1 + alpha2);
    float inv_alpha2 = 1.0f /(alpha3 + alpha4);
    float inv_alpha = 1.0f /(alpha1 + alpha2 + alpha3 + alpha4);
    float delta = 0.25 * inv_alpha1 + 0.25 * inv_alpha2;
    //gaussian product center
    float px = (alpha1 *  center1[0] + alpha2 * center2[0]) * (inv_alpha1);
    float py = (alpha1 *  center1[1] + alpha2 * center2[1]) * (inv_alpha1);
    float pz = (alpha1 *  center1[2] + alpha2 * center2[2]) * (inv_alpha1);
    float qx = (alpha3 *  center3[0] + alpha4 * center4[0]) * (inv_alpha2);
    float qy = (alpha3 *  center3[1] + alpha4 * center4[1]) * (inv_alpha2);
    float qz = (alpha3 *  center3[2] + alpha4 * center4[2]) * (inv_alpha2);
    //distances
    float pax = px -center1[0],pbx = px-center2[0];
    float pay = py -center1[1],pby = py-center2[1];
    float paz = pz -center1[2],pbz = pz-center2[2];
    float qcx = qx -center3[0],qdx = qx-center4[0];
    float qcy = qy -center3[1],qdy = qy-center4[1];
    float qcz = qz -center3[2],qdz = qz-center4[2];

    float qpdx = px-qx;
    float qpdy = py-qy;
    float qpdz = pz-qz;
    float pqd2 = qpdx * qpdx + qpdy * qpdy + qpdz * qpdz;
    //squared distance AB
    float distx1 = center1[0] -center2[0];
    float disty1 = center1[1] -center2[1];
    float distz1 = center1[2] -center2[2];
    float sq_dist1 = distx1*distx1 +disty1*disty1 +distz1*distz1;
    float distx2 = center3[0] - center4[0];
    float disty2 = center3[1] - center4[1];
    float distz2 = center3[2] - center4[2];
    float sq_dist2 = distx2*distx2 +disty2*disty2 +distz2*distz2;

    float prefact = 2 * M_PI * M_PI * inv_alpha1 * inv_alpha2 * sqrt(M_PI * inv_alpha);
    prefact *= exp(-alpha1 * alpha2 * inv_alpha1 * sq_dist1); 
    prefact *= exp(-alpha3 * alpha4 * inv_alpha2 * sq_dist2);
    float accum = 0.0f;
    for (int lp = 0; lp < l1 +l2 + 1; lp++){
    for (int rp = 0; rp < (int)(lp / 2) + 1; rp++){
    for (int lq = 0; lq < l3 +l4 + 1; lq++){
    for (int rq = 0; rq < (int)(lq / 2) + 1; rq++){
        for (int i = 0; i < (int)((lp+lq-2*rp-2*rq)/2) + 1; i++){
            float gx = gfn(lp,lq,rp,rq,i,l1,l2,l3,l4,pax,pbx,qcx,qdx,qpdx,inv_alpha1,inv_alpha2,delta);
            for (int mp = 0; mp < m1+m2 +1 ; mp++){
            for (int sp = 0; sp < (int)(mp / 2) + 1; sp++){
            for (int mq = 0; mq < m3+m4 +1 ; mq++){
            for (int sq = 0; sq < (int)(mq / 2) + 1; sq++){
                for (int j = 0; j < (int)((mp+mq-2*sp-2*sq)/2) + 1 ; j++){
                    float gy = gfn(mp,mq,sp,sq,j,m1,m2,m3,m4,pay,pby,qcy,qdy,qpdy,inv_alpha1,inv_alpha2,delta);
                    for (int np = 0; np < n1+n2 + 1; np++){
                    for (int tp = 0; tp < (int)(np / 2) + 1; tp++){
                    for (int nq = 0; nq < n3+n4 + 1; nq++){
                    for (int tq = 0; tq < (int)(nq / 2) + 1; tq++){
                        for (int k = 0; k < (int)((np+nq-2*tp-2*tq)/2) + 1; k++){
                            float gz = gfn(np,nq,tp,tq,k,n1,n2,n3,n4,paz,pbz,qcz,qdz,qpdz,inv_alpha1,inv_alpha2,delta);
                           
                            int v = lp + lq + mp + mq + np + nq 
                                - 2 *(rp + rq + sp + sq + tp + tq) - (i + j + k);
                            float F = boys(v,pqd2 / (4*delta));
                            accum += gx * gy * gz * F;
                        }
                    }}}}
                }
            }}}}
        }
    }}}}
    
    return prefact * accum;
}

void gto_electron_electron(cgto * cgtos, int cgto_num, float * G){
    for (int i = 0; i < cgto_num; i++){
        for (int j = 0; j < cgto_num; j++){
            for (int k = 0; k < cgto_num; k++){
                for (int l = 0; l < cgto_num; l++){
                    float g = 0.0f;
                    for (int ip = 0; ip < cgtos[i].n_prim; ip++){
                        for (int jp = 0; jp < cgtos[j].n_prim; jp++){
                            for (int kp = 0; kp < cgtos[k].n_prim; kp++){
                                for (int lp = 0; lp < cgtos[l].n_prim; lp++){
                                    g +=prim_electron_ee(&(cgtos[i].primitives[ip]),&(cgtos[j].primitives[jp]),
                                        &(cgtos[k].primitives[kp]),&(cgtos[l].primitives[lp]),
                                        cgtos[i].center,cgtos[j].center,cgtos[k].center,cgtos[l].center) *
                                        cgtos[i].primitives[ip].coeff * cgtos[j].primitives[jp].coeff *
                                        cgtos[k].primitives[kp].coeff * cgtos[l].primitives[lp].coeff *
                                        cgtos[i].primitives[ip].normf * cgtos[j].primitives[jp].normf * 
                                        cgtos[k].primitives[kp].normf * cgtos[l].primitives[lp].normf ;
                                   
                                }
                            }
                        }
                    }
                    G[l + cgto_num * (k + cgto_num * (j + cgto_num * i))] = g;
                }
            }
        }
    }
}