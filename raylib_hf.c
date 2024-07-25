#include "raylib_hf.h"

void hf_init(){
    return ;
};

void scf_render(void){
    return ;
};
void hf_geom_render(atoms * ats, int number){
    for (int i = 0; i < number; i++)
    {
        DrawCircleGradient(
            ats[i].pos[0] * CELL_TO_WIDHT,
            ats[i].pos[1] * CELL_TO_HEIGHT,
            (CELL_DEPTH+0.1 - ats[i].pos[2]) * 10,
            WHITE, DARKGRAY
        );
    }
    
    return ;
};
void scf_update(float dt){
    return ;
};
void hf_geom_update(float dt){
    return ;
};
