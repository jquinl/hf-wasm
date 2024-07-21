#include "hf_wasm.h"

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
            ats[i].pos.x * CELL_TO_WIDHT,
            ats[i].pos.y * CELL_TO_HEIGHT,
            (CELL_DEPTH+0.1 - ats[i].pos.z) * 10,
            WHITE, DARKGRAY
        );
    }
    
    return ;
};
void scf_update(f32 dt){
    return ;
};
void hf_geom_update(f32 dt){
    return ;
};
