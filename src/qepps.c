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
#include <math.h>

typedef struct
{
  PetscBool Active;
  Mat Matrix;
} BaseMat;

PetscErrorCode loadSweepParameters( PetscInt *nParams, PetscComplex vec_params[] )
{
  int N,p;
  float value;
  
  /* ------------------------------------------------ */
  
  PetscBool flg;
  PetscInt i;
  
  FILE *fp;
  char filename[PETSC_MAX_PATH_LEN];
  
  /* ------------------------------------------------ */
  
  PetscOptionsGetString(NULL,"-params",filename,PETSC_MAX_PATH_LEN,&flg);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must supply parameter sweep file in input arguments");
  PetscFOpen(PETSC_COMM_SELF,filename,"r",&fp);
  if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Unable to open file specifying sweep parameters");
  
  /* ------------------------------------------------ */
  
  N=0;
  while ( fscanf(fp,"%*f") != EOF )
  {
    N++;
  }
  
  PetscMalloc(N*sizeof(PetscComplex),&vec_params);
  PetscPrintf(PETSC_COMM_WORLD,"INFO: Found %d parameter values...\n",N);
  *nParams=N;
  
  rewind(fp);
  p=0;
  while( fscanf(fp,"%f",&value) != EOF )
  {
    vec_params[p]=value;
    PetscPrintf(PETSC_COMM_WORLD,"param=%f\n",vec_params[p]);
    p++;
  }
  
  PetscFClose(PETSC_COMM_SELF,fp);
  
  return 0;
}

PetscErrorCode loadMatricies( const char *optStringArray[], BaseMat baseMatrixArray[], const PetscInt baseMatrixArraySize )
{
  PetscInt i;
  PetscViewer viewer;
  PetscBool flg;
  PetscErrorCode ierr;
  
  FILE *fp;
  char filename[PETSC_MAX_PATH_LEN];
  
  for(i=0; i<=baseMatrixArraySize-1; i++ )
  {
    PetscOptionsGetString(NULL,optStringArray[i],filename,PETSC_MAX_PATH_LEN,&flg);
    if (flg) 
    {
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);
      MatCreate(PETSC_COMM_WORLD,&baseMatrixArray[i].Matrix);
      MatSetFromOptions(baseMatrixArray[i].Matrix);
      MatLoad(baseMatrixArray[i].Matrix,viewer);
      PetscViewerDestroy(&viewer);
      baseMatrixArray[i].Active=1;
    } 
    else
    {
      baseMatrixArray[i].Active=0;
    }
  }
  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  /* ------------------------------------------------ */
  
  QEP     qep;          
  QEPType type;
  
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
  
  PetscPrintf(PETSC_COMM_WORLD,"INFO: Setting up solver...\n");
  QEPCreate(PETSC_COMM_WORLD,&qep);
  QEPSetFromOptions(qep);
  
  /* ------------------------------------------------ */
  
  if( Eb[0].Active )
    MatDuplicate( Eb[0].Matrix, MAT_COPY_VALUES, &E );
  
  if( Db[0].Active )
    MatDuplicate( Db[0].Matrix, MAT_COPY_VALUES, &D );

  if( Kb[0].Active )
    MatDuplicate( Kb[0].Matrix, MAT_COPY_VALUES, &K );
    
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
    QEPSetTarget(qep,lambda_tgt);
    
    /* ------------------------------------------------ */
    
    PetscPrintf(PETSC_COMM_WORLD,"INFO: Running solver...\n");
    QEPSolve(qep);
    QEPGetConverged(qep,&nConverged);
    QEPGetIterationNumber(qep,&nIterations);
    PetscPrintf(PETSC_COMM_WORLD,"INFO: Number of iterations: %d\n",nIterations);
    
    for (ev=0; ev<nConverged; ev++)
    {
      QEPGetEigenpair( qep, ev, &lambda_solved, NULL, NULL, NULL );
      QEPComputeRelativeError( qep, ev, &error );
      PetscPrintf(PETSC_COMM_WORLD,"      lambda=%f\n",lambda_solved);
      //QEPGetEigenvector( qep, ev, Ur, Ui );
    }
    lambda_tgt=lambda_solved;
  }
  
  /* ------------------------------------------------ */
  
  PetscFree(vec_params);
  PetscFree(vec_lambdas);
  
  QEPDestroy(&qep);
  MatDestroy(&E);
  MatDestroy(&D);
  MatDestroy(&K);
  
  for (i=0; i<=2; i++)
  {
    MatDestroy(&Eb[i].Matrix);
    MatDestroy(&Db[i].Matrix);
    MatDestroy(&Kb[i].Matrix);
  }
  
  SlepcFinalize();
  
  /* ------------------------------------------------ */
  
  return 0;
}

