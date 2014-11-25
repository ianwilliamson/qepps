#include <slepcpep.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "lcomplex.h"
#include "luavars.h"
#include "types.h"

double complex getComplexNumberLUA(lua_State *L)
{
  double complex value;
  if ( lua_type(L,-1) == LUA_TNUMBER )
  {
    value=lua_tonumber(L,-1)+I*0;
  }
  else if( lua_type(L,-1) == LUA_TUSERDATA )
  {
    value=*( (double complex *)lua_touserdata(L,-1) );
  }
  return value;
}

PetscComplex getPetscComplexLUA(lua_State *L)
{
  double complex value=getComplexNumberLUA(L);
  return TO_PETSC_COMPLEX(value);
}

int getAraryLengthLUA(lua_State *L,const char* array_name)
{
  lua_getglobal(L,array_name);
  int N = lua_rawlen(L, -1);
  lua_pop(L,1);
  return N; 
}

lua_State *openConfigLUA(const char* filename_settings)
{  
  lua_State *L=NULL;
  L=luaL_newstate();
  luaL_openlibs(L);
  luaL_requiref(L, "complex", &luaopen_complex, 1);
  
  if ( luaL_dofile(L, filename_settings) )
  {
    error(L, "Error with configuration file: %s", lua_tostring(L, -1));
    exit(1);
  }
  
  return L;
}

ParameterSet *parseConfigParametersLUA(lua_State *L)
{
  int i, N;
  
  N = getAraryLengthLUA(L,LUA_var_parameters);
  ParameterSet *parameters = malloc( PARAMETER_SET_SIZE(N) );
  if (parameters==NULL)
  {
    PetscPrintf(PETSC_COMM_WORLD,"Parameter set allocation failed! \n");
    exit(1);
  }
  
  parameters->num=N;
  
  lua_getglobal(L,LUA_var_parameters);
  if (lua_istable(L, -1))
  {
    for(i=1; i<=N; i++)
    {
      lua_rawgeti(L, -1, i);
      parameters->param[i-1]=getPetscComplexLUA(L);
      lua_pop(L,1);
    }
  }
  lua_pop(L,1);
  
  return parameters;
}

MatrixComponent *parseConfigMatrixLUA(lua_State *L, const char* array_name)
{
  int i, N;
  char filename[PETSC_MAX_PATH_LEN];
  PetscViewer viewer;
  
  N = getAraryLengthLUA(L,array_name);
  MatrixComponent *M = malloc( MATRIX_COMPONENT_SIZE(N) );
  if (M==NULL)
  {
    PetscPrintf(PETSC_COMM_WORLD,"Matrix component '%s' allocation failed!\n",array_name);
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
      
      lua_pop(L,1);
    }
  }
  lua_pop(L,1);
  
  return M;
}

void deleteMatrix(MatrixComponent *M)
{
  int i;
  for(i=0; i < M->num; i++)
  {
    MatDestroy( &(M->matrix[i]) );
  }
  free(M);
}

