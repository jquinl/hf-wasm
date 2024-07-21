#include <math.h>
#include "hf_wasm.h"
//Trygve Helgaker, Poul Jorgensen, Jeppe Olsen - Molecular Electronic-Structure Theory

float gto_overlap(gto * gto1, gto * gto2){
    float inv_alpha = 1 /(gto1->alpha1 + gto2->alpha2);
    float fact = pow(M_PI * inv_alpha,1.5);
    float dx = gto1->center.x -gto1->center.x;
    float dy = gto1->center.y -gto1->center.y;
    float dz = gto1->center.z -gto1->center.z;
    
    float sq_dist = dx*dx +dy*dy +dz*dz; 
    return fact * exp(-gto1->alpha1 * gto2->alpha2 * sq_dist  / inv_alpha);
}