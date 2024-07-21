#ifndef GAME_H_
#define GAME_H_

#include <raylib.h>

#define WIDTH  600
#define HEIGHT 600
#define CELL_WIDTH  10
#define CELL_HEIGHT 10
#define CELL_DEPTH  10

#define CELL_TO_WIDHT  (WIDTH /CELL_WIDTH)
#define CELL_TO_HEIGHT (HEIGHT /CELL_HEIGHT)

#define GTOS_PER_ORBITAL 2
#define ORBITALS_PER_ATOM 2

typedef struct gto
{
    Vector3 center;
    float norm;
    float alpha;
} gto;

typedef struct orbital
{
    Vector3 center;
    gto gtos[GTOS_PER_ORBITAL];
} orbital;

typedef struct atoms
{
    int type;
    Vector3 pos;
    orbital * orbitals;//orbitals[ORBITALS_PER_ATOM];
    int num_orbitals;
} atoms;

void hf_init();

void scf_render(void);
void hf_geom_render(atoms * ats, int number);

void scf_update(f32 dt);
void hf_geom_update(f32 dt);

#endif // GAME_H_