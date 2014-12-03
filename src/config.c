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

#include <lua.h>
#include <lauxlib.h>
#include <petscmat.h>
#include <stdarg.h>
#include "types.h"
#include "luavars.h"
#include "lcomplex.h"
#include "log.h"

double complex returnComplexLUA(lua_State *L)
{
  double complex result;
  if ( lua_type(L,-1) == LUA_TNUMBER )
    result=lua_tonumber(L,-1)+I*0;
  else if( lua_type(L,-1) == LUA_TUSERDATA )
    result=*( (double complex *)lua_touserdata(L,-1) );
  else
    logError("#! LUA: Requested option is not of type 'double complex'\n");
  return result;
}

static void pullFromTableLUA(lua_State *L,const char *table,const char *option)
{
  lua_getglobal(L,table);
  if (!lua_istable(L, -1))
    logError("#! LUA: '%s' is not a table\n", lua_tostring(L, -1));
  lua_pushstring(L, option);
  lua_gettable(L, -2);
}

static void pullFromArrayLUA(lua_State *L,const char *array,int index)
{
  index++; // LUA is 1-indexed
  lua_getglobal(L,array);
  if (!lua_istable(L, -1))
    logError("#! LUA: '%s' is not an array\n", lua_tostring(L, -1));
  lua_rawgeti(L, -1, index);
}

double complex getParameterValue(lua_State *L,int index)
{
  double complex result;
  pullFromArrayLUA(L,LUA_array_parameters,index);
  result=returnComplexLUA(L);
  lua_pop(L,2);
  return result;
}

char *getOptStringLUA(lua_State *L,const char *option,const char *default_value)
{
  char *result;
  pullFromTableLUA(L,LUA_array_options,option);
  if (lua_isstring(L, -1)) {
    result=strdup(lua_tostring(L,-1));
    lua_pop(L,2); //pop string and table
  } else {
    result=strdup(default_value);
    logOutput("# LUA: '%s[%s]' is not a string, using default: %s\n",LUA_array_options,option,default_value);
    lua_pop(L,1); //pop table
  }
  return result;
}

bool getOptBooleanLUA(lua_State *L,const char *option, bool default_value)
{
  bool result;
  pullFromTableLUA(L,LUA_array_options,option);
  if (lua_isboolean(L, -1)) {
    result=lua_toboolean(L,-1);
    lua_pop(L,2); //pop boolean and table
  } else {
    result=default_value;
    logOutput("# LUA: '%s[%s]' is not a boolean, using default: %i \n",LUA_array_options,option,default_value);
    lua_pop(L,1); //pop table
  }
  return result;
}

double complex getOptComplexLUA(lua_State *L,const char *option,double complex default_value)
{
  double complex result;
  pullFromTableLUA(L,LUA_array_options,option);
  if ( lua_type(L,-1) == LUA_TNUMBER ) {
    result=lua_tonumber(L,-1)+I*0;
    lua_pop(L,2); //pop value and table
  } else if( lua_type(L,-1) == LUA_TUSERDATA ) {
    result=*( (double complex *)lua_touserdata(L,-1) );
    lua_pop(L,2); //pop value and table
  } else {
    result=default_value;
    logOutput("# LUA: '%s[%s]' is not a double complex, using default: %f%+fj\n",creal(result),cimag(result));
    lua_pop(L,1); //pop table
  }
  return result;
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
    logError("#! Error with configuration file: %s", lua_tostring(L, -1));
    /* Use basic error handling as the logger isn't setup at this point */

  return L;
}

MatrixComponent *parseConfigMatrixLUA(lua_State *L, const char* matrix_key)
{
  PetscInt m,n;
  int ne, nm, Nelems, Nmats, nf;
  char filename[PETSC_MAX_PATH_LEN];
  PetscViewer viewer;
  
  lua_getglobal(L,LUA_table_matricies);
  if ( !lua_istable(L,-1) )
    logError("#! LUA: '%s' is not a table\n",LUA_table_matricies);
  
  lua_pushstring(L,matrix_key);
  lua_gettable(L,-2);
  if ( !lua_istable(L,-1) )
    logError("#! LUA: '%s[%s]' is not a table\n",LUA_table_matricies,matrix_key);
    
  Nelems = lua_rawlen(L, -1);
  if(Nelems % 2)
    logError("#! LUA: Length of '%s[%s]' is odd\n",LUA_table_matricies,matrix_key);
  
  Nmats = Nelems/2;
  MatrixComponent *M = malloc( MATRIX_COMPONENT_SIZE(Nmats) );
  if (M==NULL)
    logError("#! Allocation of MatrixComponent container specified by '%s[%s]' failed\n",LUA_table_matricies,matrix_key);
  M->num = Nmats;
  
  nm=0; // Matrix index
  nf=0; // Function index
  ne=1; // Element index
  while(ne<=Nelems)
  {
    lua_rawgeti(L, -1, ne);
    
    if( lua_type(L,-1) == LUA_TSTRING ) {
      strcpy(filename, lua_tostring(L,-1));
      PetscViewerBinaryOpen( PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer );
      MatCreate( PETSC_COMM_WORLD, &(M->matrix[nm]) );
      MatSetType( M->matrix[nm], MATMPIAIJ );      
      MatLoad( M->matrix[nm], viewer );
      PetscViewerDestroy( &viewer );
      
      MatGetSize(M->matrix[nm],&m,&n);
      logOutput("# %dx%d matrix loaded from '%s'\n",m,n,filename);
      lua_pop(L,1);
      nm++;
      ne++;
    } else if( lua_type(L,-1) == LUA_TFUNCTION ) {
      lua_pop(L,1);
      nf++; // Count number of functions
      ne++;
    } else {
      lua_pop(L,1);
      ne++;
    }
  }
  lua_pop(L,2);

  if ( nf != nm ) // We must have a scaling function for each matrix that we've loaded
    logError("#! LUA: Mismatch in the number of functions and datafiles specified in '%s[%s]'\n",LUA_table_matricies,matrix_key);
  
  return M;
}

void deleteMatrix(MatrixComponent *M)
{
  int i;
  for(i=0; i < M->num; i++) {
    MatDestroy( &(M->matrix[i]) );
  }
  free(M);
}

