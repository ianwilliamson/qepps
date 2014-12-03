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
// Functions for outputing to stdout and a file simultaneously
// 
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <grvy.h>
#include <slepcsys.h>
#include "types.h"

static FILE *fp=NULL;

static bool isRank0(void)
{
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  return rank==0;
}

void logError(const char *format, ...)
{
  if(isRank0())
  {
    va_list args;
    va_start(args,format);
    vprintf(format,args);
    if(fp!=NULL)
    {
      va_start(args,format);
      vfprintf(fp,format,args);
      fclose(fp);
      fp=NULL;
    }
    va_end(args);
  }
  exit(1);
}

void logOutput(const char *format, ...)
{
  if(isRank0())
  {
    va_list args;
    va_start(args,format);
    vprintf(format,args);
    if(fp!=NULL)
    {
      va_start(args,format);
      vfprintf(fp,format,args);
    }
    va_end(args);
  }
}

void logOpen(const char *filename)
{
  if(isRank0() && fp==NULL)
  {
    grvy_check_file_path(filename);
    fp=fopen(filename,"w");
  }
}

void logClose(void)
{
  if(isRank0() && fp!=NULL)
  {
    fclose(fp);
    fp=NULL;
  }
}

