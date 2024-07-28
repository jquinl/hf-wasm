#ifndef PRIMITIVES_H_
#define PRIMITIVES_H_

#include <math.h>
#include "functions.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846f
#endif

typedef struct prim
{
    float normf;
    float coeff;
    float alpha;
    int l;
    int m;
    int n;
} prim;
typedef struct cgto
{
    float center[3];
    prim * primitives;
    int n_prim;
} cgto;

typedef struct atoms
{
    int type;
    float pos[3];
} atoms;


prim build_primitive(float coef,float alpha,int l, int m, int n);
int build_cgto(int at_num, atoms  * ats, cgto * cgtos,prim * primitives);
#endif