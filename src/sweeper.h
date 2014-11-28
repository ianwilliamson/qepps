#ifndef QEPPS_SWEEPER
#define QEPPS_SWEEPER

/*!
 *  This assembles the system matricies for each parameter value and handles all function evaluation
 *  within the LUA state. Handles scaling/combining the component matricies.
 */
void qeppsSweeper(lua_State *L);

/*!
 *  This is the main driver function of the QEPPS package. Assumes that the LUA state has 
 *  been primed with some configuration file.
 */
void assembleMatrix(lua_State *L, const char* array_name, Mat M, MatrixComponent *Mc, int p);

#endif
