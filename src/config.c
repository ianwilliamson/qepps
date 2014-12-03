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

static lua_State *L=NULL;

static double complex returnComplexLUA()
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

static void pullFromTableLUA(const char *table,const char *option)
{
  lua_getglobal(L,table);
  if (!lua_istable(L, -1))
    logError("#! LUA: '%s' is not a table\n", lua_tostring(L, -1));
  lua_pushstring(L, option);
  lua_gettable(L, -2);
}

static void pullFromArrayLUA(const char *array,int index)
{
  lua_getglobal(L,array);
  if (!lua_istable(L, -1))
    logError("#! LUA: '%s' is not an array\n", lua_tostring(L, -1));
  lua_rawgeti(L, -1, index+1);
}

double complex getParameterValue(int index)
{
  double complex result;
  pullFromArrayLUA(LUA_array_parameters,index);
  if ( lua_type(L,-1) == LUA_TNUMBER ) {
    result=lua_tonumber(L,-1)+I*0;
    lua_pop(L,2); //pop value and table
  } else if( lua_type(L,-1) == LUA_TUSERDATA ) {
    result=*( (double complex *)lua_touserdata(L,-1) );
    lua_pop(L,2); //pop value and table
  } else {
    result=0;
  }
  return result;
}

char *getOptStringLUA(const char *option,const char *default_value)
{
  char *result;
  pullFromTableLUA(LUA_array_options,option);
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

bool getOptBooleanLUA(const char *option, bool default_value)
{
  bool result;
  pullFromTableLUA(LUA_array_options,option);
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

double complex getOptComplexLUA(const char *option,double complex default_value)
{
  double complex result;
  pullFromTableLUA(LUA_array_options,option);
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

int getAraryLengthLUA(const char* array_name)
{
  lua_getglobal(L,array_name);
  int N = lua_rawlen(L, -1);
  lua_pop(L,1);
  return N; 
}

int getNumberOfParameters()
{
  return getAraryLengthLUA(LUA_array_parameters);
}

void startLUA(void)
{
  if(L==NULL) {
    L=luaL_newstate();
    luaL_openlibs(L);
    luaL_requiref(L, "complex", &luaopen_complex, 1);
  }
}  

void closeLUA(void)
{
  if(L!=NULL) {
    lua_close(L); 
    L=NULL;
  }
}

void parseConfigLUA(const char* filename)
{
  if ( luaL_dofile(L, filename) )
    logError("#! Error parsing config file: %s\n", lua_tostring(L, -1));
}

double complex funcParamValue(const char* matrix_name, int p, int m)
{
  double complex result=0;
  
  lua_getglobal(L,LUA_table_matricies);
  if ( !lua_istable(L,-1) )
    logError("#! LUA: '%s' is not a table\n",LUA_table_matricies);
  
  lua_pushstring(L,matrix_name);
  lua_gettable(L,-2);
  if ( !lua_istable(L,-1) )
    logError("#! LUA: '%s[%s]' is not a table\n",LUA_table_matricies,matrix_name);
  
  lua_pushstring(L,LUA_subkey_func);
  lua_gettable(L,-2);
  if ( !lua_istable(L,-1) )
    logError("#! LUA: '%s[%s][%s]' is not a table\n",LUA_table_matricies,matrix_name,LUA_subkey_func);
  
  lua_getglobal(L,LUA_array_parameters);
  if (!lua_istable(L, -1))
    logError("#! LUA: '%s' is not an array\n", lua_tostring(L, -1));
  
  lua_rawgeti(L,-2,m+1);
  if( lua_type(L,-1) == LUA_TFUNCTION) {
    lua_rawgeti(L,-2,p+1);
    lua_call(L, 1, 1); // call function
    result=returnComplexLUA(); // read value
    lua_pop(L,2); //Pop value returned value and param array
  } else {
    logError("#! LUA: Non-function type found in '%s[%s][%s]'\n",LUA_table_matricies,matrix_name,LUA_subkey_func);
  }
  lua_pop(L,2); //Pop func array and parent matrix array
  return result;
}

MatrixComponent *parseConfigMatrixLUA(const char* matrix_name)
{
  int Nfiles, Nfuncs, i, m, n;
  char filename[PETSC_MAX_PATH_LEN];
  PetscViewer viewer;
  
  lua_getglobal(L,LUA_table_matricies);
  if ( !lua_istable(L,-1) )
    logError("#! LUA: '%s' is not a table\n",LUA_table_matricies);
  
  lua_pushstring(L,matrix_name);
  lua_gettable(L,-2);
  if ( !lua_istable(L,-1) )
    logError("#! LUA: '%s[%s]' is not a table\n",LUA_table_matricies,matrix_name);
  
  lua_pushstring(L,LUA_subkey_data);
  lua_gettable(L,-2);
  if ( !lua_istable(L,-1) )
    logError("#! LUA: '%s[%s][%s]' is not a table\n",LUA_table_matricies,matrix_name,LUA_subkey_data);
  
  lua_pushstring(L,LUA_subkey_func);
  lua_gettable(L,-3);
  if ( !lua_istable(L,-1) )
    logError("#! LUA: '%s[%s][%s]' is not a table\n",LUA_table_matricies,matrix_name,LUA_subkey_func);
  
  Nfiles = lua_rawlen(L, -2);
  Nfuncs = lua_rawlen(L, -1);
  lua_pop(L,1); // Pop func table
  
  if(Nfiles!=Nfuncs)
    logError("#! LUA: Mismatch in number of functions and data files\n",LUA_table_matricies,matrix_name);
  
  MatrixComponent *M = malloc( MATRIX_COMPONENT_SIZE(Nfiles) );
  if (M==NULL)
    logError("#! Allocation of MatrixComponent container for '%s' failed\n",matrix_name);
  M->num = Nfiles;
  
  for(i=0; i < M->num; i++)
  {
    lua_rawgeti(L, -1, i+1);
    if( lua_type(L,-1) == LUA_TSTRING ) {
      strcpy(filename, lua_tostring(L,-1));
      PetscViewerBinaryOpen( PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer );
      MatCreate( PETSC_COMM_WORLD, &(M->matrix[i]) );
      MatSetType( M->matrix[i], MATMPIAIJ );      
      MatLoad( M->matrix[i], viewer );
      PetscViewerDestroy( &viewer );
      
      MatGetSize(M->matrix[i],&m,&n);
      logOutput("# %dx%d matrix loaded from '%s'\n",m,n,filename);
    } else {
      logError("#! LUA: Non-string type found in '%s[%s][%s]'\n",LUA_table_matricies,matrix_name,LUA_subkey_data);
    }
    lua_pop(L,1);
  }
  lua_pop(L,2);
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

