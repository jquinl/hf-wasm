#include "raylib_hf.h"

void hf_init(int at_num,atoms * ats, prim * prims, cgto * cgtos, int  * cgto_n,float * S,float * T,float * V,float * EE){

    *cgto_n =  build_cgto(at_num,ats,cgtos,prims);
    //Cleanup matrices
    for (int i = 0; i < *cgto_n* *cgto_n; i++){
        S[i] = 0.0f;
        T[i] = 0.0f;
        V[i] = 0.0f;
        for (int j = 0; j < *cgto_n * *cgto_n; j++){
            EE[j+i* *cgto_n * *cgto_n] = 0.0f;
        }
    }
    
    gto_overlap(cgtos, *cgto_n,S);
    gto_kinetic(cgtos, *cgto_n,T);
    gto_electron_nuclear(cgtos, *cgto_n,ats,at_num,V);
    gto_electron_electron(cgtos, *cgto_n,EE);
};

void scf_render(int cgto_n,int n_elect, float * Sinv,float * T,float * V,float * EE,float * P,float * TMP){
  //  for (int i = 0; i < cgto_n; i++){
  //      for (int j = 0; j < cgto_n; j++){
  //          DrawText(TextFormat("%f", S[i + j*10]), 100 * i, 100 * j, 10, BLACK);
  //      }
  //  }
};
void hf_geom_render(atoms * ats, int number,int selected,int hover){
    for (int i = 0; i < number; i++){
        bool hovered = hover - 1 == i ;
        Color hov = hovered ? GREEN : WHITE;
        Color col = selected-1 == i ? RED : hov;
        DrawCircleGradient(
            ats[i].pos[0] * CELL_TO_WIDTH,
            ats[i].pos[1] * CELL_TO_HEIGHT,
            (CELL_DEPTH+0.1 - ats[i].pos[2]) * 10,
            col, DARKGRAY
        );
    }
    for(int i = 0 ; i<number; i++){
        bool hovered = hover - 1 == i ;
        if (hovered){
            for (int j = 0; j < number; j++){
                if (i==j) continue;
                DrawLine(
                    ats[i].pos[0] * CELL_TO_WIDTH,
                    ats[i].pos[1] * CELL_TO_HEIGHT,
                    ats[j].pos[0] * CELL_TO_WIDTH,
                    ats[j].pos[1] * CELL_TO_HEIGHT,
                    BLACK);
                Vector2 c = {
                    0.5f * (ats[i].pos[0]+ats[j].pos[0]) * CELL_TO_WIDTH,
                    0.5f * (ats[i].pos[1]+ats[j].pos[1]) * CELL_TO_HEIGHT
                    };
                DrawCircleLinesV(c, 5.0f , BLACK);
                float d1 = ats[i].pos[0] -ats[j].pos[0];
                float d2 = ats[i].pos[1] -ats[j].pos[1];
                float d = sqrtf(d1*d1+d2*d2);
                DrawText(TextFormat("[%.2f B]", d), (int)c.x, (int)c.y + 5, 10, BLACK);
            }
        }
    }
    
    return ;
};
float scf_update( int cgto_n,int n_elect, float * Sinv,float * T,float * V,float * EE,float * P,float * TMP){

    return restricted_hartree_fock_step(Sinv,T,V,EE,cgto_n,n_elect,P,TMP);
};
void hf_geom_update(atoms * ats, int number){
    return ;
};
