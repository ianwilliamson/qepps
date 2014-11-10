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

#include "slepcqep.h"
#include "const_qepps.h"
#include "load.h"
#include <math.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  /* ------------------------------------------------ */
  
  QEP     qep;          
  QEPType type;
  ST      st;
  
  PetscMPIInt    rank;
  PetscViewer    viewer;
  PetscBool      flg;
  
  FILE *fp;
  char filename[PETSC_MAX_PATH_LEN];
  
  /* ------------------------------------------------ */
  
  PetscComplex lambda_solved, lambda_tgt, *vec_lambdas, *vec_params;
  PetscReal    error, tol;
  PetscInt     p, i, ev, nConverged, maxIterations, nIterations, nParams;
  
  /* ------------------------------------------------ */
  
  Vec Ur, Ui; // soln vectors
  
  BaseMat Eb[3], Db[3], Kb[3]; // base matricies, loaded from disk
  Mat     E,     D,     K;     // complete matricies, storing scaled values 
  
  const char *optsE[3]= {"-E0","-E1","-E2"}; // Input option strings so we can loop
  const char *optsD[3]= {"-D0","-D1","-D2"};
  const char *optsK[3]= {"-K0","-K1","-K2"};
  
  /* ------------------------------------------------ */
  
  /* INITIALIZE */
  SlepcInitialize(&argc,&argv,(char*)0,help);
  
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  /* LOAD */
  PetscPrintf(PETSC_COMM_WORLD,"INFO: Loading problem data...\n");
  
  /* FREQUENCIES */
  loadSweepParameters(&nParams,vec_params);
  
  /* MATRICIES */
  PetscPrintf(PETSC_COMM_WORLD,"INFO: Loading E...\n");
  loadMatricies( optsE, Eb, 3 );
  PetscPrintf(PETSC_COMM_WORLD,"INFO: Loading D...\n");
  loadMatricies( optsD, Db, 3 );
  PetscPrintf(PETSC_COMM_WORLD,"INFO: Loading K...\n");
  loadMatricies( optsK, Kb, 3 );
  
  /* INIT */
  PetscMalloc(nParams*sizeof(PetscComplex),&vec_lambdas);
  PetscOptionsGetScalar(NULL,"-lambda_tgt",&lambda_tgt,&flg);
  if (!flg)
    lambda_tgt=1;
  PetscPrintf(PETSC_COMM_WORLD,"INFO: lambda_tgt = %f + j%f\n",PetscRealPart(lambda_tgt),PetscImaginaryPart(lambda_tgt));
  
  /* ------------------------------------------------ */
  
  if( Eb[0].Active )
    MatDuplicate( Eb[0].Matrix, MAT_COPY_VALUES, &E );
  
  if( Db[0].Active )
    MatDuplicate( Db[0].Matrix, MAT_COPY_VALUES, &D );

  if( Kb[0].Active )
    MatDuplicate( Kb[0].Matrix, MAT_COPY_VALUES, &K );
    
  QEPCreate(PETSC_COMM_WORLD,&qep);
  QEPSetFromOptions(qep); 
  
  for (p=0; p<=nParams-1; p++)
  {
    for (i=1; i<=2; i++)
    {
      if( Eb[i].Active )
      {
        if( Eb[i-1].Active ) {
          MatAXPY( E, pow(vec_params[p],i), Eb[i].Matrix, DIFFERENT_NONZERO_PATTERN );
        } else {
          MatDuplicate( Eb[i].Matrix, MAT_COPY_VALUES, &E );
          MatScale( Eb[i].Matrix , pow(vec_params[p],i) );
        }
      }
      if( Db[i].Active )
      {
        if( Db[i-1].Active ) {
          MatAXPY( D, pow(vec_params[p],i), Db[i].Matrix, DIFFERENT_NONZERO_PATTERN );
        } else {
          MatDuplicate( Db[i].Matrix, MAT_COPY_VALUES, &D );
          MatScale( Db[i].Matrix , pow(vec_params[p],i) );
        }
      }
      if( Kb[i].Active )
      {
        if( Kb[i-1].Active ) {
          MatAXPY( K, pow(vec_params[p],i), Kb[i].Matrix, DIFFERENT_NONZERO_PATTERN );
        } else {
          MatDuplicate( Kb[i].Matrix, MAT_COPY_VALUES, &K );
          MatScale( Kb[i].Matrix , pow(vec_params[p],i) );
        }
      }
    }

    QEPSetOperators(qep,E,D,K);
    
    /* ------------------------------------------------ */
    
    PetscPrintf(PETSC_COMM_WORLD,"INFO: Running solver...\n");
    
    QEPSetTarget(qep,lambda_tgt);
    
    //QEPGetST(qep,&st);
    //STSetShift(st,lambda_tgt);
    
    QEPSolve(qep);
    QEPGetConverged(qep,&nConverged);
    QEPGetIterationNumber(qep,&nIterations);
    
    for (ev=0; ev<nConverged; ev++)
    {
      QEPGetEigenpair( qep, ev, &lambda_solved, NULL, NULL, NULL );
      QEPComputeRelativeError( qep, ev, &error );
      PetscPrintf(PETSC_COMM_WORLD,"      lambda = %f%+fj\n",PetscRealPart(lambda_solved),PetscImaginaryPart(lambda_solved));
      lambda_tgt=lambda_solved;
    }
  }
  
  /* ------------------------------------------------ */
  
  //PetscFree(vec_params);
  //PetscFree(vec_lambdas);
  
  QEPDestroy(&qep);
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

