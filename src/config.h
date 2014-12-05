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
// lua script file, and loading the data matricies, etc
// 
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef QEPPS_CONFIG
#define QEPPS_CONFIG

/*! 
 *  Returns parameters[index+1] from the LUA state
 */
double complex getParameterValue(int index);

/*! 
 *  Returns the size of the parameters table in the LUA state
 */
int getNumberOfParameters();

/*! 
 *  Returns a string from the QEPPS options table in the LUA state, returns default_value
 *  if option is not defined
 *  
 *  free() must be called on the returned pointer
 */
char *getOptStringLUA(const char *option,const char *default_value);

/*! 
 *  Returns a bolean from the QEPPS options table in the LUA state, returns default_value
 *  if option is not defined
 */
bool getOptBooleanLUA(const char *option, bool default_value);

/*! 
 *  Returns a complex double from the QEPPS options table in the LUA state, returns default_value
 *  if option is not defined
 */
double complex getOptComplexLUA(const char *option,double complex default_value);

/*! 
 *  Returns an int from the QEPPS options table in the LUA state, returns default_value
 *  if option is not defined
 */
int getOptIntLUA(const char *option,int default_value);

/*! 
 *  Returns the length of the LUA array identified by the string array_name. Pushes and pops
 *  from the stack so the stack should be in the same state as before the call.
 */
int getAraryLengthLUA(const char* array_name);

/*!
 *  Sets up a new LUA state and opens the default LUA librariesas well as the complex number
 *  library.
 */
void startLUA(void);

/*!
 *  Closes the LUA state.
 */
void closeLUA(void);

/*!
 *  Runs the configuration script identified by filename.
 */
void parseConfigLUA(const char* filename);

/*!
 *  Returns the result of evaluating the m-th function of the 'matrix_name' component on the p-th
 *  parameter value. 
 *  
 *  p and m are indexed from zero (C style rather than LUA style)
 */
double complex funcParamValue(const char* matrix_name, int p, int m);

/*!
 *  Parses and loads the matricies from the data files specified in LUA
 */
MatrixComponent *parseConfigMatrixLUA(const char* matrix_name);

/*!
 *  Traverses the MatrixComponent struct and calls MatDestroy on each of the listed Mat's. After
 *  this, calls free() on the MatrixComponent struct.
 */
void deleteMatrix(MatrixComponent *M);

#endif
