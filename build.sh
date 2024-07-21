#!/bin/bash
emcc -Os -Wall  -I. -I${RAYLIB_PATH} -L. -L${RAYLIB_PATH} -l:libraylib.a -s USE_GLFW=3 -s ASYNCIFY --shell-file shell.html -DPLATFORM_WEB -o out/hf_wasm.html *.c 