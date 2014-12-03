//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QEPPS: Quadratic eigenvalue problem parameter sweeper
//
// Copyright (C) 2014 Lab for Active Nano Devices, UT ECE 
// Developed by Ian Williamson 
// Supervised by Dr. Zheng Wang 
//
//-----------------------------------------------------------------------el-
// 
// Configuration subsystem for interacting with the LUA state, loading the
// lua script file, and loading the data matricies
// 
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef QEPPS_CONFIG
#define QEPPS_CONFIG

/*! 
 *  
 */
double complex getParameterValue(lua_State *L,int index);

/*! 
 *  Returns a string from the QEPPS options table in the LUA state
 *  Note that free() must be called
 */
char *getOptStringLUA(lua_State *L,const char *option,const char *default_value);

/*! 
 *  Returns a bolean from the QEPPS options table in the LUA state
 */
bool getOptBooleanLUA(lua_State *L,const char *option, bool default_value);

/*! 
 *  Returns a complex double from the QEPPS options table in the LUA state
 */
double complex getOptComplexLUA(lua_State *L,const char *option,double complex default_value);

/*!
 *  Returns the most recently pushed variable on the LUA stack as a complex double data type.
 *  This assumes that the value to be returned has already been pushed onto the stack, 
 *  i.e. with a call to lua_getglobal() or as the result of function eval
 *  Note: this function does not pop the value from the stack.
 */
double complex returnComplexLUA(lua_State *L);

/*! 
 *  Returns the length of the LUA array identified by the string array_name. Pushes and pops
 *  from the stack so the stack should be in the same state as before the call.
 */
int getAraryLengthLUA(lua_State *L,const char* array_name);

/*!
 *  Sets up a new LUA state and opens the default LUA libraries Also opens the complex numbers
 *  library. Runs the configuration script identified by the string filename_settings.
 *  Returns the LUA state.
 */
lua_State *openConfigLUA(const char* filename_settings);

/*!
 *  Parses, from the specified LUA state, the matrix data specified in the array identified
 *  by the input string array_name.
 */
MatrixComponent *parseConfigMatrixLUA(lua_State *L, const char* array_name);

/*!
 *  Traverses the MatrixComponent struct and calls MatDestroy on each of the listed Mat's. After
 *  this it frees the entire MatrixComponent struct.
 */
void deleteMatrix(MatrixComponent *M);

#define getNumberOfParameters(L) getAraryLengthLUA(L,LUA_array_parameters)

#endif
