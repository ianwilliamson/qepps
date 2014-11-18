#include "slepcpep.h"
#include "const_qepps.h"

PetscErrorCode loadSweepParameters( PetscInt *nParams, PetscReal **vec_params )
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
    N++;
  
  PetscMalloc( ((size_t)N)*sizeof(PetscReal), &vec_params );
  *nParams=N;
  
  rewind(fp);
  p=0;
  while( fscanf(fp,"%f",&value) != EOF )
  {
    (*vec_params)[p]=(PetscReal)value;
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
    baseMatrixArray[i].Active=0;
    
    PetscOptionsGetString(NULL,optStringArray[i],filename,PETSC_MAX_PATH_LEN,&flg);
    if (flg) 
    {
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);
      MatCreate(PETSC_COMM_WORLD,&baseMatrixArray[i].Matrix);
      MatSetType(baseMatrixArray[i].Matrix,MATMPIAIJ);      
      MatLoad(baseMatrixArray[i].Matrix,viewer);
      PetscViewerDestroy(&viewer);
      
      baseMatrixArray[i].Active=1;
    }
  }
  return 0;
}
