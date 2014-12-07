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

static void assembleMatrix(const char* matrix_name, Mat M, MatrixComponent *Mc, int p)
{
  int i;
  double complex value;
  
  MatZeroEntries(M);
  
  for(i=0;i<Mc->num;i++)
  {
    value=funcParamValue(matrix_name,p,i);
    MatAXPY( M, TO_PETSC_COMPLEX(value), Mc->matrix[i], DIFFERENT_NONZERO_PATTERN );
  }
  
  MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);  
}

void qeppsSweeper(void)
{
  PEP pep;  
  ST st;     
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
  MatrixComponent *Ec = parseConfigMatrixLUA(LUA_key_matrix_E);
  MatrixComponent *Dc = parseConfigMatrixLUA(LUA_key_matrix_D);
  MatrixComponent *Kc = parseConfigMatrixLUA(LUA_key_matrix_K);
  
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
  lambda_tgt = getOptComplexLUA("lambda_tgt",1);
  logOutput("# lambda_tgt set to %.3f%+.3fj\n",creal(lambda_tgt),cimag(lambda_tgt));
  
  // Initialize the solver
  A[0]=K; A[1]=D; A[2]=E;
  PEPCreate(PETSC_COMM_WORLD,&pep);
  PEPSetProblemType(pep,PEP_GENERAL);
  PEPGetST(pep,&st);
  STSetTransform(st,1);
  STSetType(st,STSINVERT);
  PEPSetDimensions(pep,getOptIntLUA("nev",1),PETSC_DEFAULT,PETSC_DEFAULT);
  PEPSetFromOptions(pep);
  
  MPI_Comm_size(PETSC_COMM_WORLD,&p); 
  logOutput("# MPI_Comm_size = %i \n", p);
  logOutput("# Number of parameters = %i \n", getNumberOfParameters());
  grvy_timer_end("setup");
  for (p=0; p < getNumberOfParameters(); p++)
  {
    grvy_timer_begin("assemble");
    logOutput("%E", getParameterValue(p) );
    
    assembleMatrix(LUA_key_matrix_E,E,Ec,p);
    assembleMatrix(LUA_key_matrix_D,D,Dc,p);
    assembleMatrix(LUA_key_matrix_K,K,Kc,p);
    
    MatGetVecs(E,&Uout,NULL);
    MatGetVecs(E,&Uinit,NULL);
    
    PEPSetOperators(pep,3,A);
    PEPSetTarget(pep,TO_PETSC_COMPLEX(lambda_tgt));
    grvy_timer_end("assemble");
    
    grvy_timer_begin("solve");
    PEPSolve(pep);
    grvy_timer_end("solve");
    
    grvy_timer_begin("postprocess");
    PEPGetConverged(pep,&nConverged);
    for (ev=0; ev<nConverged; ev++)
    {
      PEPGetEigenpair( pep, ev, &lambda_solved, NULL, Uout, NULL );
      logOutput(", %.3f%+.3fj",PetscRealPart(lambda_solved),PetscImaginaryPart(lambda_solved));
      
      if(ev==0) // Leading eigenvalue/eigenvector (should be closest to target)
      {
        if( getOptBooleanLUA("update_lambda_tgt", false) )
        {
          lambda_tgt = TO_DOUBLE_COMPLEX(lambda_solved);
        }
        if( getOptBooleanLUA("update_initspace", false) )
        {
          VecCopy(Uout,Uinit);
          PEPSetInitialSpace(pep,1,&Uinit);
        }
      }
      if( getOptBooleanLUA("save_solutions", false) )
      {
        char filename[PETSC_MAX_PATH_LEN];
        char *output_dir = getOptStringLUA("output_dir","./");
        sprintf(filename,"%s/U_%E_%i.dat",output_dir,creal( getParameterValue(p) ),ev);
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
    grvy_timer_end("postprocess");
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
  
  if( getOptBooleanLUA("print_timing",false) )
  {
    logOutput("# ================================================\n");
    logOutput("# ================================================\n");
    logOutput("# total time: %10.5E secs\n",grvy_timer_elapsed_global());
    logOutput("# ------------------------------------------------\n");
    logOutput("#      setup: %10.5E secs\n",grvy_timer_elapsedseconds("setup"));
    logOutput("#   assemble: %10.5E secs\n",grvy_timer_elapsedseconds("assemble"));
    logOutput("#      solve: %10.5E secs\n",grvy_timer_elapsedseconds("solve"));
    logOutput("#   postproc: %10.5E secs\n",grvy_timer_elapsedseconds("postprocess"));
    logOutput("#      clean: %10.5E secs\n",grvy_timer_elapsedseconds("clean"));
    logOutput("# ------------------------------------------------\n");    
    logOutput("# assemble  (   count): %i\n",     grvy_timer_stats_count("assemble"));
    logOutput("# assemble  (    mean): %E secs\n",grvy_timer_stats_mean("assemble"));
    logOutput("# assemble  (variance): %E secs\n",grvy_timer_stats_variance("assemble"));
    logOutput("# assemble  (     min): %E secs\n",grvy_timer_stats_min("assemble"));
    logOutput("# assemble  (     max): %E secs\n",grvy_timer_stats_max("assemble"));
    logOutput("# solve     (   count): %i\n",     grvy_timer_stats_count("solve"));
    logOutput("# solve     (    mean): %E secs\n",grvy_timer_stats_mean("solve"));
    logOutput("# solve     (variance): %E secs\n",grvy_timer_stats_variance("solve"));
    logOutput("# solve     (     min): %E secs\n",grvy_timer_stats_min("solve"));
    logOutput("# solve     (     max): %E secs\n",grvy_timer_stats_max("solve"));
    logOutput("# postproc  (   count): %i\n",     grvy_timer_stats_count("postprocess"));
    logOutput("# postproc  (    mean): %E secs\n",grvy_timer_stats_mean("postprocess"));
    logOutput("# postproc  (variance): %E secs\n",grvy_timer_stats_variance("postprocess"));
    logOutput("# postproc  (     min): %E secs\n",grvy_timer_stats_min("postprocess"));
    logOutput("# postproc  (     max): %E secs\n",grvy_timer_stats_max("postprocess"));    
    logOutput("# ================================================\n");
    logOutput("# ================================================\n");
  }
}
