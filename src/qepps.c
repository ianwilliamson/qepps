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

#include "slepcpep.h"
#include "const_qepps.h"
#include "load.h"
#include "assemble.h"
#include <math.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  /* ------------------------------------------------ */
  
  PEP     pep;          
  PEPType type;
  
  PetscMPIInt    rank;
  PetscViewer    viewer;
  PetscBool      flg;
  
  FILE *fp;
  char filename[PETSC_MAX_PATH_LEN];
  
  /* ------------------------------------------------ */
  
  PetscComplex lambda_solved, lambda_tgt, *vec_lambdas;
  PetscReal    *vec_params;
  PetscReal    error, tol;
  PetscInt     p, i, ev, nConverged, maxIterations, nIterations, nParams;
  
  /* ------------------------------------------------ */
  
  Vec Ur, Ui; // soln vectors
  
  BaseMat Eb[3], Db[3], Kb[3];           // base matricies, loaded from disk
  Mat     E,     D,     K,     A[3];     // complete matricies, storing scaled values 
  
  const char *optsE[3]= {"-E0","-E1","-E2"}; // Input option strings so we can loop
  const char *optsD[3]= {"-D0","-D1","-D2"};
  const char *optsK[3]= {"-K0","-K1","-K2"};
  
  /* ------------------------------------------------ */
  
  /* INITIALIZE */
  SlepcInitialize(&argc,&argv,(char*)0,help);
  
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  /* FREQUENCIES */
  //loadSweepParameters(&nParams,&vec_params);
  nParams=1;
  /* MATRICIES */
  loadMatricies( optsE, Eb, 3 );
  loadMatricies( optsD, Db, 3 );
  loadMatricies( optsK, Kb, 3 );
  
  /* INIT */
  //PetscMalloc(nParams*sizeof(PetscComplex),&vec_lambdas);
  PetscOptionsGetScalar(NULL,"-lambda_tgt",&lambda_tgt,&flg);
  if (!flg) lambda_tgt=1;
  
  /* ------------------------------------------------ */
  MatCreate(PETSC_COMM_WORLD,&E);
  MatCreate(PETSC_COMM_WORLD,&D);
  MatCreate(PETSC_COMM_WORLD,&K);
  MatSetType(E,MATMPIAIJ);
  MatSetType(D,MATMPIAIJ);
  MatSetType(K,MATMPIAIJ);
      
  if( Eb[0].Active )
    MatDuplicate( Eb[0].Matrix, MAT_COPY_VALUES, &E );
  
  if( Db[0].Active )
    MatDuplicate( Db[0].Matrix, MAT_COPY_VALUES, &D );

  if( Kb[0].Active )
    MatDuplicate( Kb[0].Matrix, MAT_COPY_VALUES, &K );
  
  PEPCreate(PETSC_COMM_WORLD,&pep);
  
  for (p=0; p<=nParams-1; p++)
  {
    //PetscPrintf(PETSC_COMM_WORLD,"%.3f,  ",vec_params[p]);
    for (i=1; i<=2; i++)
    {
      //incorporateMatrixComponent( E, pow( vec_params[p], i), Eb[i], Eb[i-1] );
      //incorporateMatrixComponent( D, pow( vec_params[p], i), Db[i], Db[i-1] );
      //incorporateMatrixComponent( K, pow( vec_params[p], i), Kb[i], Kb[i-1] );
    }
     
    A[0]=K; A[1]=D; A[2]=E;
    PEPSetOperators(pep,3,A);
    PEPSetProblemType(pep,PEP_GENERAL);
    PEPSetFromOptions(pep);
    
    /* ------------------------------------------------ */
    
    PEPSetTarget(pep,lambda_tgt);
    
    PEPSolve(pep);
    PEPGetConverged(pep,&nConverged);
    
    if(nConverged >= 1)
    {
      for (ev=0; ev<nConverged; ev++)
      {
        PEPGetEigenpair( pep, ev, &lambda_solved, NULL, NULL, NULL );
        PetscPrintf(PETSC_COMM_WORLD,"%.3f%+.3fj,  ",PetscRealPart(lambda_solved),PetscImaginaryPart(lambda_solved));
        if(ev==0)
          lambda_tgt=lambda_solved; // Set target for next parameter value
      }
      PetscPrintf(PETSC_COMM_WORLD,"\n");
    }
    else
    {
      // No eigen values found ...
    }
  }
  
  /* ------------------------------------------------ */
  
  //PetscFree(vec_params);
  //PetscFree(vec_lambdas);
  
  PEPDestroy(&pep);
  MatDestroy(&E);
  MatDestroy(&D);
  MatDestroy(&K);
  
  for (i=0; i<=2; i++)
  {
    if( Eb[i].Active )
      MatDestroy(&(Eb[i].Matrix));
    if( Db[i].Active )
      MatDestroy(&(Db[i].Matrix));
    if( Kb[i].Active )
      MatDestroy(&(Kb[i].Matrix));
  }
  
  SlepcFinalize();
  
  /* ------------------------------------------------ */
  
  return 0;
}

