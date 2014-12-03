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
  
  /* Start LUA state and load configuration */
  startLUA();
  PetscOptionsGetString(NULL,"-lua",filename,PETSC_MAX_PATH_LEN,NULL);
  parseConfigLUA(filename);
  
  /* Setup the log file */
  char *output_file = getOptStringLUA("output_log","./output.txt");
  logOpen(output_file);
  free(output_file);
  
  /* Run parameter sweep */
  qeppsSweeper();
  
  /* Close out */
  closeLUA();
  logClose();
  SlepcFinalize();
   
  return 0;
}
