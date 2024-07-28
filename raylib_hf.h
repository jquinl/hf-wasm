
#ifndef RAYLIB_HF_H_
#define RAYLIB_HF_H_

#include <raylib.h>
#include "hf/primitives.h"
#include "hf/integrals_naive.h"

#define WIDTH  600
#define HEIGHT 600
#define CELL_WIDTH  10
#define CELL_HEIGHT 10
#define CELL_DEPTH  10

#define CELL_TO_WIDHT  (WIDTH /CELL_WIDTH)
#define CELL_TO_HEIGHT (HEIGHT /CELL_HEIGHT)

void hf_init(int at_num,atoms * ats, prim * prims, cgto * cgtos,int * cgto_n,float * S,float * T,float * V,float * EE);

void scf_render(int at_num,atoms * ats, cgto * cgtos,int  cgto_n,float * S,float * T,float * V,float * EE);
void hf_geom_render(atoms * ats, int number);

void scf_update(int at_num,atoms * ats, cgto * cgtos,int  cgto_n,float * S,float * T,float * V,float * EE);
void hf_geom_update(atoms * ats, int number);
#endif