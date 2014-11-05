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

PetscErrorCode loadSweepParameters( PetscInt *nParams, PetscReal *vec_params )
{
  PetscBool flg;
  PetscInt i;
  PetscErrorCode ierr;
  
  FILE *fp;
  char filename[PETSC_MAX_PATH_LEN];

  ierr = PetscOptionsGetString(NULL,"-params",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must supply parameter sweep file in input arguments");

  PetscFOpen(PETSC_COMM_SELF,filename,"r",&fp);
  if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Unable to open file specifying sweep parameters");

  while ( EOF != fscanf(fp,"%*E") ) //scan but discard with *
    ++nParams;
  PetscMalloc(2*sizeof(PetscReal),&vec_params);
  PetscPrintf(PETSC_COMM_WORLD,"Found %d parameter values...\n",nParams);

  rewind(fp);
  i=0;
  while( EOF != fscanf(fp,"%E",&vec_params[i]) )
    ++i;
  PetscFClose(PETSC_COMM_SELF,fp);
  
  return 0;
}

PetscErrorCode loadMatricies( char *optStringArray[], Mat *matrixArray[], PetscInt numMatricies )
{
  PetscInt i;
  PetscViewer viewer;
  PetscBool flg;
  PetscErrorCode ierr;
  
  FILE *fp;
  char filename[PETSC_MAX_PATH_LEN];
  
  for(i=0; i<=numMatricies-1; i++ )
  {
    ierr = PetscOptionsGetString(NULL,optStringArray[i],filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
    if (flg) 
    {
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = MatCreate(PETSC_COMM_WORLD,&matrixArray[i]);CHKERRQ(ierr);
      ierr = MatSetFromOptions(matrixArray[i]);CHKERRQ(ierr);
      ierr = MatLoad(matrixArray[i],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    } 
    else
    {
      matrixArray[i]=NULL;
      //SETERRQ(PETSC_COMM_WORLD,1,"Must indicate a file name for matrix M with the -M option");
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
  PetscErrorCode ierr;
  
  FILE *fp;
  char filename[PETSC_MAX_PATH_LEN];
  
  /* ------------------------------------------------ */
  
  PetscComplex lambda_solved, lambda_tgt, *vec_lambdas;
  PetscScalar  error, tol, *vec_params, lambda_tgt_real=1, lambda_tgt_imag=0;
  PetscInt     i, j, nev, nconv, maxit, its, nParams;
  
  /* ------------------------------------------------ */
  
  Vec Ur, Ui; // soln vectors
  
  Mat *E_base[3], *D_base[3], *K_base[3];   // base matricies, loaded from disk
  Mat  E,          D,          K;           // complete matricies, storing scaled values 
  
  const char *optsE[3]= {"-E0","-E1","-E2"}; // Input option strings so we can loop
  const char *optsD[3]= {"-D0","-D1","-D2"};
  const char *optsK[3]= {"-K0","-K1","-K2"};
  
  /* ------------------------------------------------ */
  
  /* INITIALIZE */
  SlepcInitialize(&argc,&argv,(char*)0,help);
  
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  
  /* LOAD */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Loading problem data...\n");CHKERRQ(ierr);
  
  /* FREQUENCIES */
  ierr = loadSweepParameters(&nParams,&vec_params);CHKERRQ(ierr);
  
  /* MATRICIES */
  ierr = loadMatricies( optsE, E_base, 3 );CHKERRQ(ierr);
  ierr = loadMatricies( optsD, D_base, 3 );CHKERRQ(ierr);
  ierr = loadMatricies( optsK, K_base, 3 );CHKERRQ(ierr);
  
  /* SETUP */
  // allocate storage
  PetscMalloc(nParams*sizeof(PetscReal),&vec_lambdas);
  // get starting eigenvalue target
  PetscOptionsGetReal(NULL,"-lambda_tgt_real",&lambda_tgt_real,NULL);
  PetscOptionsGetReal(NULL,"-lambda_tgt_imag",&lambda_tgt_imag,NULL);
  lambda_tgt = lambda_tgt_real+PETSC_i*lambda_tgt_imag;
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting up solver...\n");CHKERRQ(ierr);
  ierr = QEPCreate(PETSC_COMM_WORLD,&qep);CHKERRQ(ierr);
  ierr = QEPSetFromOptions(qep);CHKERRQ(ierr);
  
  /* ------------------------------------------------ */
  
  for (i=0; i<=nParams-1; i++)
  {
    ierr = MatDuplicate( E_base[0], MAT_COPY_VALUES, &E );
    ierr = MatAXPY( E, vec_params[i], E_base[1], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
    ierr = MatAXPY( E, vec_params[i]*vec_params[i], E_base[2], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
    
    ierr = MatDuplicate( D_base[0], MAT_COPY_VALUES, &D );
    ierr = MatAXPY( D, vec_params[i], D_base[1], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
    ierr = MatAXPY( D, vec_params[i]*vec_params[i], D_base[2], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
    
    ierr = MatDuplicate( K_base[0], MAT_COPY_VALUES, &K );
    ierr = MatAXPY( K, vec_params[i], K_base[1], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
    ierr = MatAXPY( K, vec_params[i]*vec_params[i], K_base[2], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);

    ierr = QEPSetOperators(qep,E,D,K);CHKERRQ(ierr);
    ierr = QEPSetTarget(qep,lambda_tgt);CHKERRQ(ierr);
    
    /* ------------------------------------------------ */
    
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Running solver...\n");CHKERRQ(ierr);
    
    ierr = QEPSolve(qep);CHKERRQ(ierr);
    ierr = QEPGetConverged(qep,&nconv);CHKERRQ(ierr);
    ierr = QEPGetIterationNumber(qep,&its);CHKERRQ(ierr);
    
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Num iterations: %D\n",its);CHKERRQ(ierr);
    for (j=0; j<nconv; j++)
    {
      ierr = QEPGetEigenpair( qep, j, &lambda_tgt, NULL, NULL, NULL );CHKERRQ(ierr);
      ierr = QEPComputeRelativeError( qep, j, &error );CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"lambda=:%D\n",its);CHKERRQ(ierr);
      //ierr = QEPGetEigenvector( qep, j, Ur, Ui );CHKERRQ(ierr);
    }
  }
  
  /* ------------------------------------------------ */
  
  PetscFree(vec_params);
  PetscFree(vec_lambdas);
  
  ierr = QEPDestroy(&qep);CHKERRQ(ierr);
  ierr = MatDestroy(&E);CHKERRQ(ierr);
  ierr = MatDestroy(&D);CHKERRQ(ierr);
  ierr = MatDestroy(&K);CHKERRQ(ierr);
  
  ierr = MatDestroy(&E_base[0]);CHKERRQ(ierr);
  ierr = MatDestroy(&E_base[1]);CHKERRQ(ierr);
  ierr = MatDestroy(&E_base[2]);CHKERRQ(ierr);
  
  ierr = MatDestroy(&D_base[0]);CHKERRQ(ierr);
  ierr = MatDestroy(&D_base[1]);CHKERRQ(ierr);
  ierr = MatDestroy(&D_base[2]);CHKERRQ(ierr);
  
  ierr = MatDestroy(&K_base[0]);CHKERRQ(ierr);
  ierr = MatDestroy(&K_base[1]);CHKERRQ(ierr);
  ierr = MatDestroy(&K_base[2]);CHKERRQ(ierr);
  
  ierr = SlepcFinalize();
  
  /* ------------------------------------------------ */
  
  return 0;
}

