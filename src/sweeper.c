#include <slepcpep.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include "types.h"
#include "luavars.h"
#include "config.h"

static void assembleMatrix(lua_State *L, const char* matrix_key, Mat M, MatrixComponent *Mc, int p)
{
  int Nelems;
  p++; // LUA arrays are indexed from 1
  double complex value;
  
  MatZeroEntries(M);
  
  lua_getglobal(L,LUA_array_parameters);
  if ( !lua_istable(L,-1) )
    error(L,"#! LUA: '%s' is not a table\n",LUA_array_parameters);
  
  lua_getglobal(L,LUA_table_matricies);
  if ( !lua_istable(L,-1) )
    error(L,"#! LUA: '%s' is not a table\n",LUA_table_matricies);
  
  lua_pushstring(L,matrix_key);
  lua_gettable(L,-2);
  if ( !lua_istable(L,-1) )
    error(L,"#! LUA: '%s[%s]' is not a table\n",LUA_table_matricies,matrix_key);
  
  Nelems = lua_rawlen(L, -1);
  if(Nelems % 2)
    error(L,"#! LUA: Length of '%s[%s]' is odd\n",LUA_table_matricies,matrix_key);
  
  int nm=0; // Function index
  int ne=1; // Element index
  while(ne<=Nelems && nm < Mc->num)
  {
    lua_rawgeti(L, -1, ne);
    
    if( lua_type(L,-1) == LUA_TFUNCTION) {
      lua_rawgeti(L, -4, p); // get parameter value
      lua_call(L, 1, 1); // call function
        
      value=returnComplexLUA(L); // read value
      lua_pop(L,1); // remove read value
        
      MatAXPY( M, TO_PETSC_COMPLEX(value), Mc->matrix[nm], DIFFERENT_NONZERO_PATTERN );
        
      nm++;
      ne++;
    } else {
      lua_pop(L,1);
      ne++;
    }
  }
  lua_pop(L,3); //pop all arrays 
  
  MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);  
}

static void saveSolutionVector(lua_State *L, Vec U, int p, int ev)
{
  return;
}

void qeppsSweeper(lua_State *L)
{
  PEP pep;       
  Vec U;
  Mat E, D, K, A[3];
  PetscComplex lambda_solved;
  PetscReal    error, tol;
  PetscInt     i, ev, nConverged, maxIterations, nIterations;
  int p;
  double complex lambda_tgt;
  
  // From LUA state, get and load matrix components
  MatrixComponent *Ec = parseConfigMatrixLUA(L,LUA_key_matrix_E);
  MatrixComponent *Dc = parseConfigMatrixLUA(L,LUA_key_matrix_D);
  MatrixComponent *Kc = parseConfigMatrixLUA(L,LUA_key_matrix_K);
  
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
  lambda_tgt = getOptComplexLUA(L,"lambda_tgt");
  PetscPrintf(PETSC_COMM_WORLD,"# lambda_tgt set to %.3f%+.3fj\n",creal(lambda_tgt),cimag(lambda_tgt));
  
  // Initialize the solver
  A[0]=K; A[1]=D; A[2]=E;
  PEPCreate(PETSC_COMM_WORLD,&pep);
  PEPSetProblemType(pep,PEP_GENERAL);
  PEPSetFromOptions(pep);
  
  PetscPrintf(PETSC_COMM_WORLD,"# Sweeping %d parameters\n", getNumberOfParameters(L));
  for (p=0; p < getNumberOfParameters(L); p++)
  {
    PetscPrintf(PETSC_COMM_WORLD,"%E", getParameterValue(L,p) );
    
    assembleMatrix(L,LUA_key_matrix_E,E,Ec,p);
    assembleMatrix(L,LUA_key_matrix_D,D,Dc,p);
    assembleMatrix(L,LUA_key_matrix_K,K,Kc,p);
    
    PEPSetOperators(pep,3,A);
    PEPSetTarget(pep,TO_PETSC_COMPLEX(lambda_tgt));
    PEPSolve(pep);
    PEPGetConverged(pep,&nConverged);
    
    if(nConverged >= 1)
    {
      for (ev=0; ev<nConverged; ev++)
      {
        if( getOptBooleanLUA(L,"save_solutions") || getOptBooleanLUA(L,"update_initspace") )
          PEPGetEigenpair( pep, ev, &lambda_solved, NULL, U, NULL );
        else
          PEPGetEigenpair( pep, ev, &lambda_solved, NULL, NULL, NULL );
        
        PetscPrintf(PETSC_COMM_WORLD,", %.3f%+.3fj",PetscRealPart(lambda_solved),PetscImaginaryPart(lambda_solved));
        
        if(ev==0)
        {
          if( getOptBooleanLUA(L,"update_lambda_tgt") )
            lambda_tgt = TO_DOUBLE_COMPLEX(lambda_solved);
          if( getOptBooleanLUA(L,"update_initspace") )
            PEPSetInitialSpace(pep,1,&U);
          if( getOptBooleanLUA(L,"save_solutions") )
            saveSolutionVector(L,U,p,ev);
        }
      } // loop converged
    }
    else
    {
      PetscPrintf(PETSC_COMM_WORLD,"\n");
      break; // Stop sweeping if we don't solve for any eigenvalues.
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");
  } // loop parameters
  
  PEPDestroy(&pep);
  MatDestroy(&E);
  MatDestroy(&D);
  MatDestroy(&K);
  deleteMatrix(Ec);
  deleteMatrix(Dc);
  deleteMatrix(Kc);
}
