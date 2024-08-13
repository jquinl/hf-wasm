#ifndef RHF_H_
#define RHF_H_

#include "functions.h"
#include "integrals_naive.h"
#include "utils.h"

void restricted_hartree_fock_init(const float * S,const float * T, const float * V,const float * EE,const int n_basis,const int n_elec,float * P,float* Sinvsqr,float * TMP);
float restricted_hartree_fock_step(const float * S,const float * T, const float * V,const float * EE,const int n_basis,const int n_elec,float * P,float * TMP);
float nuclear_repulsion(atoms * ats, int at_num);
#endif