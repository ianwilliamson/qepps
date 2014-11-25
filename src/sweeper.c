#include <slepcpep.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include "types.h"
#include "luavars.h"
#include "config.h"

/*!
 *  This assembles the system matricies for each parameter value and handles all function evaluation
 *  within the LUA state. Handles scaling/combining the component matricies.
 */
void assembleMatrix(lua_State *L, const char* array_name, Mat M, MatrixComponent *Mc, int p)
{
  int i; p++; // LUA arrays are indexed from 1
  double complex value;
  
  MatZeroEntries(M);
  
  lua_getglobal(L,LUA_var_parameters); // stack 1: param value array
  lua_getglobal(L,array_name);   // stack 2: function pointer array
  if (lua_istable(L, -1) && lua_istable(L, -2))
  {
    for(i=1; i<=Mc->num; i++)
    {
      lua_rawgeti(L, -1, i); //get current function
      if (lua_isfunction(L,-1))
      {
        lua_rawgeti(L, -3, p); //get parameter value
        lua_call(L, 1, 1);    //call function
      } // if not a function, it will be treated as a value
      value=getComplexNumberLUA(L); //read
      lua_pop(L,1); //remove read value
      MatAXPY( M, TO_PETSC_COMPLEX(value), Mc->matrix[i-1], DIFFERENT_NONZERO_PATTERN );
    }
  }
  lua_pop(L,2); // pop both arrays
  
  MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);  
}

/*!
 *  This is the main driver function of the QEPPS package. Assumes that the LUA state has 
 *  been primed with some configuration file.
 */
void qeppsSweeper(lua_State *L)
{
  PEP pep;       
  Vec Ur, Ui;
  Mat E, D, K, A[3];
  PetscComplex lambda_solved;
  PetscReal    error, tol;
  PetscInt     i, ev, nConverged, maxIterations, nIterations;
  int p;
  
  // From LUA state, get parameters that are to be swept
  ParameterSet *parameters = parseConfigParametersLUA(L);
  // From LUA state, get and load matrix components
  MatrixComponent *Ec = parseConfigMatrixLUA(L, LUA_array_Edat);
  MatrixComponent *Dc = parseConfigMatrixLUA(L, LUA_array_Ddat);
  MatrixComponent *Kc = parseConfigMatrixLUA(L, LUA_array_Kdat);
  
  // Initialize total matricies
  // (we scale/sum the component matricies from the previous step into these)
  MatCreate(PETSC_COMM_WORLD,&E);
  MatCreate(PETSC_COMM_WORLD,&D);
  MatCreate(PETSC_COMM_WORLD,&K);
  MatSetType(E,MATMPIAIJ);
  MatSetType(D,MATMPIAIJ);
  MatSetType(K,MATMPIAIJ);
  
  // Prime the total matricies with the problem size and nonzero pattern
  MatDuplicate(Ec->matrix[0],MAT_SHARE_NONZERO_PATTERN,&E);
  MatDuplicate(Dc->matrix[0],MAT_SHARE_NONZERO_PATTERN,&D);
  MatDuplicate(Kc->matrix[0],MAT_SHARE_NONZERO_PATTERN,&K);
  
  // Get the target eigenvalue from the LUA state
  lua_getglobal(L,LUA_var_lambda_tgt);
  PetscComplex lambda_tgt = getPetscComplexLUA(L);
  lua_pop(L,1);
  
  // Initialize the solver
  A[0]=K; A[1]=D; A[2]=E;
  PEPCreate(PETSC_COMM_WORLD,&pep);
  PEPSetProblemType(pep,PEP_GENERAL);
  PEPSetFromOptions(pep);
  
  for (p=0; p<parameters->num; p++)
  {
    PetscPrintf(PETSC_COMM_WORLD,"%E,  ",parameters->param[p]);
    
    assembleMatrix(L,LUA_array_Efuncs,E,Ec,p);
    assembleMatrix(L,LUA_array_Dfuncs,D,Dc,p);
    assembleMatrix(L,LUA_array_Kfuncs,K,Kc,p);
    
    PEPSetOperators(pep,3,A);
    //PEPSetInitialSpace(pep,PetscInt n,Vec *is)
    PEPSetTarget(pep,lambda_tgt);
    PEPSolve(pep);
    PEPGetConverged(pep,&nConverged);
    
    if(nConverged >= 1)
    {
      for (ev=0; ev<nConverged; ev++)
      {
        PEPGetEigenpair( pep, ev, &lambda_solved, NULL, NULL, NULL );
        PetscPrintf(PETSC_COMM_WORLD,"%.3f%+.3fj,  ",PetscRealPart(lambda_solved),PetscImaginaryPart(lambda_solved));
        //if(ev==0)
        //  lambda_tgt=lambda_solved; // Set target for next parameter value
      }
    }
    else
    {
      // No eigen values found ...
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");
  }
  
  PEPDestroy(&pep);
  MatDestroy(&E);
  MatDestroy(&D);
  MatDestroy(&K);
  deleteMatrix(Ec);
  deleteMatrix(Dc);
  deleteMatrix(Kc);
  free(parameters);
}
