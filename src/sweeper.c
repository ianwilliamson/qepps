#include "slepcpep.h"
#include <petscbag.h>
#include "common.h"

void qeppsSweeper(BaseMat Eb[], BaseMat Db[], BaseMat Kb[], PetscComplex lambda_tgt, PetscBag *bag)
{
  PEP       pep;          
  PEPType   type;
  
  Vec Ur, Ui;        // soln vectors 
  Mat E, D, K, A[3]; // complete matricies, storing scaled values
  SweepSet *sweep;

  PetscComplex lambda_solved, *vec_lambdas;
  PetscReal    error, tol;
  PetscInt     p, i, ev, nConverged, maxIterations, nIterations;

  PetscBagGetData(*bag,(void**)&sweep);
  
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
  
  for (p=0; p<sweep->size; p++)
  {
    PetscPrintf(PETSC_COMM_WORLD,"%.3f,  ",sweep->param[p]);
    for (i=1; i<=2; i++)
    {
      incorporateMatrixComponent( E, PetscPowScalar( sweep->param[p],i), Eb[i], Eb[i-1] );
      incorporateMatrixComponent( D, PetscPowScalar( sweep->param[p],i), Db[i], Db[i-1] );
      incorporateMatrixComponent( K, PetscPowScalar( sweep->param[p],i), Kb[i], Kb[i-1] );
    }
    
    A[0]=K; A[1]=D; A[2]=E;
    PEPSetOperators(pep,3,A);
    PEPSetProblemType(pep,PEP_GENERAL);
    PEPSetFromOptions(pep);
    
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
    
  PEPDestroy(&pep);
  MatDestroy(&E);
  MatDestroy(&D);
  MatDestroy(&K);

}
