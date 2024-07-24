#ifndef HF_SETTINGS_H_
#define HF_SETTINGS_H_

#ifndef M_PI
    #define M_PI 3.14159265358979323846f
#endif

#define PRIMITIVES_PER_ORBITAL 2
#define ORBITALS_PER_ATOM 2

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
    int center[3];
    prim primitives[PRIMITIVES_PER_ORBITAL];
} cgto;

typedef struct atoms
{
    int type;
    int pos[3];
    cgto orbitals[ORBITALS_PER_ATOM];
    int num_orbitals;
} atoms;



#endif 