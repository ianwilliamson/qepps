#ifndef QEPPS_SWEEPER
#define QEPPS_SWEEPER

void qeppsSweeper(lua_State *L);
void assembleMatrix(lua_State *L, const char* array_name, Mat M, MatrixComponent *Mc, int p);

#endif
