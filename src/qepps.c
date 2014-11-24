/****************************************************************************
 * Copyright (C) 2014 by Ian Williamson                                     *
 *                                                                          *
 * Quadratic Eigenvalue Problem Parameter-sweeper                           *
 *                                                                          * 
 ****************************************************************************/

#include <slepcpep.h>

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>

#include "help.h"
#include "types.h"
#include "sweeper.h"
#include "config.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  char filename[PETSC_MAX_PATH_LEN];
  SlepcInitialize(&argc,&argv,(char*)0,help);
  
  PetscOptionsGetString(NULL,"-lua",filename,PETSC_MAX_PATH_LEN,NULL);
  lua_State *L = openConfigLUA(filename);
  
  qeppsSweeper(L);
  
  lua_close(L);  
  SlepcFinalize();
  return 0;
}

