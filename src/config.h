#ifndef QEPPS_CONFIG
#define QEPPS_CONFIG

double complex getComplexNumberLUA(lua_State *L);
PetscComplex getPetscComplexLUA(lua_State *L);
void processMatrixLUA(lua_State *L,const char* Mfuncs,const char* Mdat);
int getAraryLengthLUA(lua_State *L,const char* array_name);

lua_State *openConfigLUA(const char* filename_settings);

ParameterSet *parseConfigParametersLUA(lua_State *L);
MatrixComponent *parseConfigMatrixLUA(lua_State *L, const char* array_name);
void deleteMatrix(MatrixComponent *M);

#endif
