/****************************************************************************
 * Copyright (C) 2014 by Ian Williamson                                     *
 *                                                                          *
 * Quadratic Eigenvalue Problem Parameter-sweeper                           *
 *                                                                          * 
 * Loads matrices representing quadratic eigenvalue problem and solves      *
 * (lambda^2*E+lambda*D+K)*U=0, where E, D, and K are matrices, U is a      * 
 * vector, and lambda is an eigenvalue. The matrix inputs (E, D, and K) are *
 * each be specified in components that can be combined in the form         *
 * E=E0+p*E1+p^2*E2, where p is some problem parameter (NOT the             * 
 * eigenvalue).                                                             *
 *                                                                          *
 * p represents some parameter that should be swept (typically frequency    *
 * or the free-space wavevector k_0 in electromagnetics) and is specified   *
 * as an ASCII file with one p value per line.                              *
 *                                                                          *
 * The code assembles the total E, D, and K matrices for each value of p    *
 * and computes the eigenvector (typically mode's field pattern) and        *
 * eigenvalue (typically the mode's complex effective index) closest to     *
 * some specified target.                                                   *
 *                                                                          * 
 ****************************************************************************/

#include <slepcpep.h>
#include <petscbag.h>
#include "help.h"
#include "common.h"
#include "load.h"
#include "assemble.h"
#include "sweeper.h"


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscBool flg;
  PetscReal vec_params[3]={1.,2.,3.};
  PetscComplex lambda_tgt;
  PetscInt i, nParams;
  PetscBag bag;
  
  BaseMat Eb[3], Db[3], Kb[3];           // base matricies, loaded from disk 
  
  const char *optsE[3]= {"-E0","-E1","-E2"}; // Input option strings so we can loop
  const char *optsD[3]= {"-D0","-D1","-D2"};
  const char *optsK[3]= {"-K0","-K1","-K2"};
  
  SlepcInitialize(&argc,&argv,(char*)0,help);
  
  loadSweepParameters(&bag);
  loadMatricies(optsE,Eb,3);
  loadMatricies(optsD,Db,3);
  loadMatricies(optsK,Kb,3);
  PetscOptionsGetScalar(NULL,"-lambda_tgt",&lambda_tgt,&flg);
  if (!flg) lambda_tgt=1;
  
  /* ------------------------------------------------ */
  
  qeppsSweeper(Eb,Db,Kb,lambda_tgt,&bag);
  
  /* ------------------------------------------------ */
    
  for (i=0; i<=2; i++)
  {
    if( Eb[i].Active )
      MatDestroy(&(Eb[i].Matrix));
    if( Db[i].Active )
      MatDestroy(&(Db[i].Matrix));
    if( Kb[i].Active )
      MatDestroy(&(Kb[i].Matrix));
  }
  
  PetscBagDestroy(&bag);
  
  SlepcFinalize();
  
  return 0;
}

