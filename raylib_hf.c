#include "raylib_hf.h"

void hf_init(int at_num,atoms * ats, prim * prims, cgto * cgtos, int  * cgto_n,float * S,float * T,float * V,float * EE){
    *cgto_n =  build_cgto(at_num,ats,cgtos,prims);
    gto_overlap(cgtos, *cgto_n,S);
    gto_kinetic(cgtos, *cgto_n,T);
    gto_electron_nuclear(cgtos, *cgto_n,ats,at_num,V);
    gto_electron_electron(cgtos, *cgto_n,EE);
};

void scf_render(int at_num,atoms * ats, cgto * cgtos,int  cgto_n,float * S,float * T,float * V,float * EE){
    for (int i = 0; i < cgto_n; i++){
        for (int j = 0; j < cgto_n; j++){
            DrawText(TextFormat("%f", S[i + j*10]), 100 * i, 100 * j, 10, BLACK);
        }
    }
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
void scf_update(int at_num,atoms * ats, cgto * cgtos,int  cgto_n, float * S,float * T,float * V,float * EE){
    return ;
};
void hf_geom_update(atoms * ats, int number){
    return ;
};
