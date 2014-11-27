#include <slepcpep.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include <stdlib.h>
#include "types.h"
#include "lcomplex.h"
#include "luavars.h"
#include "config.h"

/*! 
 *  
 */
void pullFromTableLUA(lua_State *L,const char *table,const char *option)
{
  lua_getglobal(L,table);
  if (!lua_istable(L, -1)) error(L, "Error is not a table: %s", lua_tostring(L, -1));
  lua_pushstring(L, option);
  lua_gettable(L, -2);
}

/*! 
 *  
 */
void pullFromArrayLUA(lua_State *L,const char *array,int index)
{
  lua_getglobal(L,array);
  if (!lua_istable(L, -1)) error(L, "Error is not an array: %s", lua_tostring(L, -1));
  lua_rawgeti(L, -1, index+1);
}

/*! 
 *  Returns a complex double from the QEPPS options table in the LUA state
 */
double complex getParameterValue(lua_State *L,int index)
{
  double complex result;
  pullFromArrayLUA(L,LUA_var_parameters,index);
  result=returnComplexLUA(L);
  lua_pop(L,2);
  return result;
}

/*! 
 *  Returns a string from the QEPPS options table in the LUA state
 *  Note that free() must be called
 */
char *getOptStringLUA(lua_State *L,const char *option)
{
  char *result;
  pullFromTableLUA(L,LUA_var_options_table,option);
  if (!lua_isstring(L, -1)) error(L,"Requested option is not a string");
  result=strdup( lua_tostring(L,-1) );
  lua_pop(L,2);
  return result;
}

/*! 
 *  Returns a bolean from the QEPPS options table in the LUA state
 */
bool getOptBooleanLUA(lua_State *L,const char *option)
{
  bool result;
  pullFromTableLUA(L,LUA_var_options_table,option);
  if (!lua_isboolean(L, -1)) error(L,"Requested option is not a boolean");
  result=lua_toboolean(L,-1);
  lua_pop(L,2);
  return result;
}

/*! 
 *  Returns a complex double from the QEPPS options table in the LUA state
 */
double complex getOptComplexLUA(lua_State *L,const char *option)
{
  double complex result;
  pullFromTableLUA(L,LUA_var_options_table,option);
  result=returnComplexLUA(L);
  lua_pop(L,2);
  return result;
}

/*!
 *  Returns the most recently pushed variable on the LUA stack as a complex double data type.
 *  This assumes that the value to be returned has already been pushed onto the stack, 
 *  i.e. with a call to lua_getglobal() or as the result of function eval
 *  Note: this function does not pop the value from the stack.
 */
double complex returnComplexLUA(lua_State *L)
{
  double complex result;
  if ( lua_type(L,-1) == LUA_TNUMBER )
    result=lua_tonumber(L,-1)+I*0;
  else if( lua_type(L,-1) == LUA_TUSERDATA )
    result=*( (double complex *)lua_touserdata(L,-1) );
  else
    error(L,"Requested option is not a complex double");
  return result;
}

/*! 
 *  Returns the length of the LUA array identified by the string array_name. Pushes and pops
 *  from the stack so the stack should be in the same state as before the call.
 */
int getAraryLengthLUA(lua_State *L,const char* array_name)
{
  lua_getglobal(L,array_name);
  int N = lua_rawlen(L, -1);
  lua_pop(L,1);
  return N; 
}

/*!
 *  Sets up a new LUA state and opens the default LUA libraries Also opens the complex numbers
 *  library. Runs the configuration script identified by the string filename_settings.
 *  Returns the LUA state.
 */
lua_State *openConfigLUA(const char* filename_settings)
{  
  lua_State *L=NULL;
  L=luaL_newstate();
  luaL_openlibs(L);
  luaL_requiref(L, "complex", &luaopen_complex, 1);
  //luaopen_complex(L);
  //lua_rawseti(L,LUA_REGISTRYINDEX,LUA_RIDX_GLOBALS);
  //lua_settop(L,0);
  //lua_pushnil(L);
  
  if ( luaL_dofile(L, filename_settings) )
  {
    error(L, "Error with configuration file: %s", lua_tostring(L, -1));
    exit(1);
  }
  
  PetscPrintf(PETSC_COMM_WORLD,"# Parsed configuration from from '%s'\n",filename_settings);
  
  return L;
}

/*!
 *  Parses, from the specified LUA state, the matrix data specified in the array identified
 *  by the input string array_name.
 */
MatrixComponent *parseConfigMatrixLUA(lua_State *L, const char* array_name)
{
  PetscInt m,n;
  int i, N;
  char filename[PETSC_MAX_PATH_LEN];
  PetscViewer viewer;
  
  N = getAraryLengthLUA(L,array_name);
  MatrixComponent *M = malloc( MATRIX_COMPONENT_SIZE(N) );
  if (M==NULL)
  {
    PetscPrintf(PETSC_COMM_WORLD,"#! Matrix component '%s' allocation failed!\n",array_name);
    exit(1);
  }
  M->num = N;
  
  lua_getglobal(L,array_name);
  if (lua_istable(L, -1))
  {
    for(i=1; i<=N; i++)
    {
      lua_rawgeti(L, -1, i);
      strcpy(filename, lua_tostring(L,-1));
      PetscViewerBinaryOpen( PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer );
      MatCreate( PETSC_COMM_WORLD, &(M->matrix[i-1]) );
      MatSetType( M->matrix[i-1], MATMPIAIJ );      
      MatLoad( M->matrix[i-1], viewer );
      PetscViewerDestroy( &viewer );
      
      MatGetSize(M->matrix[i-1],&m,&n);
      PetscPrintf(PETSC_COMM_WORLD,"# %dx%d matrix loaded from '%s'\n",m,n,filename);
      
      lua_pop(L,1);
    }
  }
  lua_pop(L,1);
  
  return M;
}

/*!
 *  Traverses the MatrixComponent struct and calls MatDestroy on each of the listed Mat's. After
 *  this it frees the entire MatrixComponent struct.
 */
void deleteMatrix(MatrixComponent *M)
{
  int i;
  for(i=0; i < M->num; i++)
  {
    MatDestroy( &(M->matrix[i]) );
  }
  free(M);
}

