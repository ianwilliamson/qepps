#ifndef QEPPS_SWEEPER
#define QEPPS_SWEEPER

/*!
 *  This assembles the system matricies for each parameter value and handles all function evaluation
 *  within the LUA state. Handles scaling/combining the component matricies.
 */
void qeppsSweeper(lua_State *L);

#endif
