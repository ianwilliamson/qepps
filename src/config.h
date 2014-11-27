#ifndef QEPPS_CONFIG
#define QEPPS_CONFIG

void pullFromTableLUA(lua_State *L,const char *table,const char *option);
void pullFromArrayLUA(lua_State *L,const char *array,int index);
double complex getParameterValue(lua_State *L,int index);
char *getOptStringLUA(lua_State *L,const char *option);
bool getOptBooleanLUA(lua_State *L,const char *option);
double complex getOptComplexLUA(lua_State *L,const char *option);
double complex returnComplexLUA(lua_State *L);
int getAraryLengthLUA(lua_State *L,const char* array_name);
lua_State *openConfigLUA(const char* filename_settings);
MatrixComponent *parseConfigMatrixLUA(lua_State *L, const char* array_name);
void deleteMatrix(MatrixComponent *M);

#define getNumberOfParameters(L) getAraryLengthLUA(L,LUA_var_parameters)

#endif
