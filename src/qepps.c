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
// Initializes Petsc and LUA and then calls the sweeper driver function
// 
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <lua.h>
#include <slepcsys.h> 
#include "help.h"
#include "types.h"
#include "sweeper.h"
#include "config.h"
#include "log.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  char filename[PETSC_MAX_PATH_LEN];
  SlepcInitialize(&argc,&argv,(char*)0,help);
  
  /* Load configuration */
  PetscOptionsGetString(NULL,"-lua",filename,PETSC_MAX_PATH_LEN,NULL);
  lua_State *L = openConfigLUA(filename);
  
  /* Setup the log file */
  char *output_file = getOptStringLUA(L,"output_log","./output.txt");
  logOpen(output_file);
  free(output_file);
  
  /* Run parameter sweep */
  qeppsSweeper(L);
  
  /* Close out */
  lua_close(L);
  logClose();
  SlepcFinalize();
   
  return 0;
}
