
#ifndef RAYLIB_HF_H_
#define RAYLIB_HF_H_

#include <raylib.h>
#include "hf/primitives.h"
#include "hf/integrals_naive.h"
#include "hf/functions.h"
#include "hf/rhf.h"

#define WIDTH  600
#define HEIGHT 600
#define CELL_WIDTH  5
#define CELL_HEIGHT 5
#define CELL_DEPTH  5

#define CELL_TO_WIDTH  (WIDTH /CELL_WIDTH)
#define CELL_TO_HEIGHT (HEIGHT /CELL_HEIGHT)

void hf_init(int at_num,atoms * ats, prim * prims, cgto * cgtos,int * cgto_n,float * S,float * T,float * V,float * EE);

void scf_render(int cgto_n,int n_elect, float * Sinv,float * T,float * V,float * EE,float * P,float * TMP);
void hf_geom_render(atoms * ats, int number,int selected,int hover);

float scf_update(int cgto_n,int n_elect, float * Sinv,float * T,float * V,float * EE,float * P,float * TMP);
void hf_geom_update(atoms * ats, int number);
#endif