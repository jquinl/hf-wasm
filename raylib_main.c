#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

#include "raylib_hf.h"

#define PLATFORM_WEB 1
#if defined(PLATFORM_WEB)
    #include <emscripten/emscripten.h>
#endif

#define MAX_NUM_CGTO 25
#define MAX_NUM_PRIM MAX_NUM_CGTO * 3
#define EPSILON 1e-5

prim prims[MAX_NUM_PRIM];
cgto cgtos[MAX_NUM_CGTO];

float S[MAX_NUM_CGTO * MAX_NUM_CGTO];
float T[MAX_NUM_CGTO * MAX_NUM_CGTO];
float V[MAX_NUM_CGTO * MAX_NUM_CGTO];
float EE[MAX_NUM_CGTO * MAX_NUM_CGTO * MAX_NUM_CGTO * MAX_NUM_CGTO];

float Sisqr[MAX_NUM_CGTO * MAX_NUM_CGTO];
float P[MAX_NUM_CGTO * MAX_NUM_CGTO];
float TMP[MAX_NUM_CGTO * MAX_NUM_CGTO * MAX_NUM_CGTO * MAX_NUM_CGTO];

//placeholder
atoms atoms_list[MAX_NUM_CGTO / 5] ;
int cgto_num = 0;
int atom_num = 0;
int n_elec = 0;
float e_nuc = 0.0f;
float last_energy = 1000.0f;
float new_energy = 1000.0f;
float delta_energy = 1000.0f;
//UI----------------------
const Color UI_RUNNING = { 193, 225, 193, 255 };
bool runSCF = false;
bool modifAt = true;
bool conv = false;
int moved_atom = 0;
int hover_atom = 0;

Rectangle run_scf_button = { WIDTH * 0.1, HEIGHT * 0.9, 100.0f, 30.0f };

void UpdateDrawFrame(void);    
int atom_hover(atoms* ats,int at_n, Vector2 mousepos);
int atom_hover_select(atoms* ats,int at_n, Vector2 mousepos);

int main(void)
{
    atom_num = 2;
    n_elec = 2;
    cgto_num = 2;
    atoms_list[0]=(atoms){1,{1.0f,1.0f,0.0f}};
    atoms_list[1]=(atoms){1,{3.0f,3.0f,0.0f}};
    
    InitWindow(WIDTH, HEIGHT, "Hartree-Fock");


#if defined(PLATFORM_WEB)
    emscripten_set_main_loop(UpdateDrawFrame, 0, 1);
#else
    SetTargetFPS(60);   
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        UpdateDrawFrame();
    }
#endif
    
    CloseWindow();

    return 0;
}

void UpdateDrawFrame(void)
{
    Vector2 mouse = GetMousePosition();

    bool is_hover = false;
    if (CheckCollisionPointRec(mouse, run_scf_button)){
        is_hover = true;
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT))
        {
            if(!runSCF && modifAt){
                hf_init(atom_num,atoms_list,prims,cgtos,&cgto_num,S,T,V,EE);
                e_nuc = nuclear_repulsion(atoms_list,atom_num);
                restricted_hartree_fock_init(S,T,V,EE,cgto_num,n_elec,P,Sisqr,TMP);
                last_energy = 1000.0f;
                conv = false;
            }
            runSCF = !runSCF;
        }
    }
    //Drag and Drop atoms
    hover_atom = atom_hover(atoms_list,atom_num,mouse);
    moved_atom = atom_hover_select(atoms_list,atom_num,mouse);
    if(!IsCursorOnScreen()) moved_atom =0;

    if(moved_atom){
        runSCF = false;
        conv = false;
        modifAt = true;
    } 
    if (moved_atom){
        atoms_list[moved_atom-1].pos[0] = (float)mouse.x / (float)CELL_TO_WIDTH;
        atoms_list[moved_atom-1].pos[1] = (float)mouse.y / (float)CELL_TO_HEIGHT;
        atoms_list[moved_atom-1].pos[2] = 0.0f;
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) moved_atom= 0;
    }
    if(runSCF){
        new_energy =e_nuc + scf_update(cgto_num,n_elec,Sisqr,T,V,EE,P,TMP);
        delta_energy = new_energy-last_energy;
        conv =  fabs(delta_energy) < 10e-3f;
        last_energy = new_energy;
        if(conv){
            runSCF = false;
            modifAt = false;
        }
    }

    BeginDrawing();
        ClearBackground(DARKGRAY);
        hf_geom_render(atoms_list,atom_num,moved_atom,hover_atom);

        DrawRectangleRec(run_scf_button,runSCF ? UI_RUNNING: LIGHTGRAY);
        DrawRectangleLines(run_scf_button.x,run_scf_button.y,run_scf_button.width,run_scf_button.height,BLACK);
        int y_ofset = is_hover ? 3 : 5;
        if(runSCF){
            if(conv){
                ///Not necessary to draw here
            }else{
                DrawText(TextFormat("Hartree-Fock Running- dE : %f",delta_energy), 10, 20, 20, BLACK);
            }
            //scf_render(2,atoms_list,cgtos,cgto_num,S,T,V,EE);
            DrawText("Stop", run_scf_button.x + run_scf_button.width/2 - MeasureText("Stop", 10), (int)run_scf_button.y +y_ofset, 20, BLACK);

        }else{
            DrawText("Run", run_scf_button.x + run_scf_button.width/2 - MeasureText("Run", 10), (int)run_scf_button.y +y_ofset, 20, BLACK);
        }
        if(conv && !runSCF){
            DrawText(TextFormat("Hartree-Fock Converged! - E : %f hartrees",last_energy), 10, 20, 20, BLACK);
        }else if(!runSCF){
            DrawText("Hartree-Fock", 10, 20, 20, BLACK);
        }
    EndDrawing();

}
int atom_hover_select(atoms* ats,int at_n, Vector2 mousepos){
    for (int i = 0; i < at_n; i++){
        Vector2 p = {ats[i].pos[0]*CELL_TO_WIDTH,ats[i].pos[1]* CELL_TO_HEIGHT};
        if (CheckCollisionPointCircle(mousepos,p, (CELL_DEPTH+0.1 - ats[i].pos[2]) * 10) && IsMouseButtonDown(MOUSE_BUTTON_LEFT)){
            return i+1;
        }
    }
    return 0;
}
int atom_hover(atoms* ats,int at_n, Vector2 mousepos){
    for (int i = 0; i < at_n; i++){
        Vector2 p = {ats[i].pos[0]*CELL_TO_WIDTH,ats[i].pos[1]* CELL_TO_HEIGHT};
        if (CheckCollisionPointCircle(mousepos,p, (CELL_DEPTH+0.1 - ats[i].pos[2]) * 10) ){
            return i+1;
        }
    }
    return 0;
}