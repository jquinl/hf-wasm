#include <raylib.h>

#include "hf_settings.h"

#define WIDTH  600
#define HEIGHT 600
#define CELL_WIDTH  10
#define CELL_HEIGHT 10
#define CELL_DEPTH  10

#define CELL_TO_WIDHT  (WIDTH /CELL_WIDTH)
#define CELL_TO_HEIGHT (HEIGHT /CELL_HEIGHT)

void hf_init();

void scf_render(void);
void hf_geom_render(atoms * ats, int number);

void scf_update(float dt);
void hf_geom_update(float dt);