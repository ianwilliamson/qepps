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
// The driver function and its auxiliaries
// 
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <lua.h>
#include <slepcpep.h>
#include <grvy.h>
#include "types.h"
#include "luavars.h"
#include "config.h"
#include "log.h"

static void assembleMatrix(lua_State *L, const char* matrix_key, Mat M, MatrixComponent *Mc, int p)
{
  int Nelems;
  p++; // LUA arrays are indexed from 1
  double complex value;
  
  MatZeroEntries(M);
  
  lua_getglobal(L,LUA_array_parameters);
  if ( !lua_istable(L,-1) )
    logError("#! LUA: '%s' is not a table\n",LUA_array_parameters);
  
  lua_getglobal(L,LUA_table_matricies);
  if ( !lua_istable(L,-1) )
    logError("#! LUA: '%s' is not a table\n",LUA_table_matricies);
  
  lua_pushstring(L,matrix_key);
  lua_gettable(L,-2);
  if ( !lua_istable(L,-1) )
    logError("#! LUA: '%s[%s]' is not a table\n",LUA_table_matricies,matrix_key);
  
  Nelems = lua_rawlen(L, -1);
  if(Nelems % 2)
    logError("#! LUA: Length of '%s[%s]' is odd\n",LUA_table_matricies,matrix_key);
  
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

static void saveSolutionVector(lua_State *L, Vec *U, int p, int ev)
{
  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vector.dat",FILE_MODE_WRITE,&viewer);
  VecView(*U,viewer);
  PetscViewerDestroy(&viewer);
  VecDestroy(U);
  return;
}

void qeppsSweeper(lua_State *L)
{
  PEP pep;       
  Vec Uout, Uinit;
  Mat E, D, K, A[3];
  PetscComplex lambda_solved;
  PetscReal    error, tol;
  PetscInt     i, ev, nConverged, maxIterations, nIterations;
  PetscViewer  viewer;
  int p;
  double complex lambda_tgt;
  
  grvy_timer_init("qepps_parameter_sweep");
  grvy_timer_begin("setup");
  
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
  lambda_tgt = getOptComplexLUA(L,"lambda_tgt",1);
  logOutput("# lambda_tgt set to %.3f%+.3fj\n",creal(lambda_tgt),cimag(lambda_tgt));
  
  // Initialize the solver
  A[0]=K; A[1]=D; A[2]=E;
  PEPCreate(PETSC_COMM_WORLD,&pep);
  PEPSetProblemType(pep,PEP_GENERAL);
  PEPSetFromOptions(pep);
  
  logOutput("# Sweeping %d parameters\n", getNumberOfParameters(L));
  grvy_timer_end("setup");
  for (p=0; p < getNumberOfParameters(L); p++)
  {
    grvy_timer_begin("iteration");
    
    logOutput("%E", getParameterValue(L,p) );
    
    assembleMatrix(L,LUA_key_matrix_E,E,Ec,p);
    assembleMatrix(L,LUA_key_matrix_D,D,Dc,p);
    assembleMatrix(L,LUA_key_matrix_K,K,Kc,p);
    
    MatGetVecs(E,&Uout,NULL);
    MatGetVecs(E,&Uinit,NULL);
    
    PEPSetOperators(pep,3,A);
    PEPSetTarget(pep,TO_PETSC_COMPLEX(lambda_tgt));
    PEPSolve(pep);
    PEPGetConverged(pep,&nConverged);
    
    for (ev=0; ev<nConverged; ev++)
    {
      PEPGetEigenpair( pep, ev, &lambda_solved, NULL, Uout, NULL );
      logOutput(", %.3f%+.3fj",PetscRealPart(lambda_solved),PetscImaginaryPart(lambda_solved));
      
      if(ev==0) // Leading eigenvalue/eigenvector (should be closest to target)
      {
        if( getOptBooleanLUA(L,"update_lambda_tgt", false) )
        {
          lambda_tgt = TO_DOUBLE_COMPLEX(lambda_solved);
        }
        if( getOptBooleanLUA(L,"update_initspace", false) )
        {
          VecCopy(Uout,Uinit);
          PEPSetInitialSpace(pep,1,&Uinit);
        }
      }
      if( getOptBooleanLUA(L,"save_solutions", false) )
      {
        char filename[PETSC_MAX_PATH_LEN];
        char *output_dir = getOptStringLUA(L,"output_dir","./");
        sprintf(filename,"%s/U_%E_%i.dat",output_dir,creal( getParameterValue(L,p) ),ev);
        free(output_dir);
        grvy_check_file_path(filename);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);
        VecView(Uout,viewer);
        PetscViewerDestroy(&viewer);
      }
    }
    logOutput("\n");
    
    if(nConverged==0)
      logError("#! Solver did not converge. Aborting...\n");
    
    grvy_timer_end("iteration");
  } // loop parameters
  
  grvy_timer_begin("clean");
  PEPDestroy(&pep);
  MatDestroy(&E);
  MatDestroy(&D);
  MatDestroy(&K);
  deleteMatrix(Ec);
  deleteMatrix(Dc);
  deleteMatrix(Kc);
  VecDestroy(&Uout);
  VecDestroy(&Uinit);
  grvy_timer_end("clean");
  
  grvy_timer_finalize();
  
  if( getOptBooleanLUA(L,"print_timing",false) )
  {
    logOutput("# \n");
    logOutput("# total time: %10.5E secs\n",grvy_timer_elapsed_global());
    logOutput("# \n");
    logOutput("#      setup: %10.5E secs\n",grvy_timer_elapsedseconds("setup"));
    logOutput("#  iteration: %10.5E secs\n",grvy_timer_elapsedseconds("iteration"));
    logOutput("#      clean: %10.5E secs\n",grvy_timer_elapsedseconds("clean"));
    logOutput("# \n");
    logOutput("# iteration (   count): %i\n",grvy_timer_stats_count("iteration"));
    logOutput("# iteration (    mean): %E secs\n",grvy_timer_stats_mean("iteration"));
    logOutput("# iteration (variance): %E secs\n",grvy_timer_stats_variance("iteration"));
    logOutput("# iteration (     min): %E secs\n",grvy_timer_stats_min("iteration"));
    logOutput("# iteration (     max): %E secs\n",grvy_timer_stats_max("iteration"));
  }
}
