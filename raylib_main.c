#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>


#include "raylib_hf.h"

#define PLATFORM_WEB 1
#if defined(PLATFORM_WEB)
    #include <emscripten/emscripten.h>
#endif

//placeholder
atoms atoms_list[2] = {
    {1,{1.0,1.0,1.0}},
    {2,{9.0,9.0,9.0}}
};

//UI----------------------
const Color UI_RUNNING = { 193, 225, 193, 255 };
bool runSCF = false;
Rectangle run_scf_button = { WIDTH * 0.1, HEIGHT * 0.9, 100.0f, 30.0f };

void UpdateDrawFrame(void);    


int main(void)
{
    InitWindow(WIDTH, HEIGHT, "Hartree-Fock");
    hf_init();

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
    bool is_hover = false;
    if (CheckCollisionPointRec(GetMousePosition(), run_scf_button)){
        is_hover = true;
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT))
        {
            runSCF = !runSCF;
        }
    }
    if(runSCF){
        
    }

    BeginDrawing();
        ClearBackground(DARKGRAY);
        hf_geom_render(atoms_list,2);
        DrawText("Hartree-Fock", 10, 20, 20, BLACK);
        DrawRectangleRec(run_scf_button,runSCF ? UI_RUNNING: LIGHTGRAY);
        DrawRectangleLines(run_scf_button.x,run_scf_button.y,run_scf_button.width,run_scf_button.height,BLACK);
        int y_ofset = is_hover ? 3 : 5;
        if(runSCF){
            DrawText("Stop", run_scf_button.x + run_scf_button.width/2 - MeasureText("Stop", 10), (int)run_scf_button.y +y_ofset, 20, BLACK);
        }else{
            DrawText("Run", run_scf_button.x + run_scf_button.width/2 - MeasureText("Run", 10), (int)run_scf_button.y +y_ofset, 20, BLACK);
        }
    EndDrawing();
}