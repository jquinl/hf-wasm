#ifndef INTEGRALS_NAIVE_H_
#define INTEGRALS_NAIVE_H_


#include "functions.h"
#include "primitives.h"

void gto_overlap(cgto * cgtos, int cgto_num, float * S);
void gto_kinetic(cgto * cgtos, int cgto_num,float * T);
void gto_electron_nuclear(cgto * cgtos, int cgto_num, atoms * atom_array, int atom_num,float * V);
void gto_electron_electron(cgto * cgtos, int cgto_num, float * G);
#endif