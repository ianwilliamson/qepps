#include "slepcpep.h"
#include <petscbag.h>
#include "common.h"

void loadSweepParameters( PetscBag *bag )
{
  int i;
  float value;  
  PetscBool flg;
  FILE *fp=NULL;
  char filename[PETSC_MAX_PATH_LEN];
  SweepSet *sweep;
  
  PetscOptionsGetString(NULL,"-params",filename,PETSC_MAX_PATH_LEN,&flg);
  if(flg) PetscFOpen(PETSC_COMM_SELF,filename,"r",&fp);
  if(fp)
  {
    for(i=0; fscanf(fp,"%*f") != EOF; i++ )
    
    PetscBagCreate(PETSC_COMM_WORLD,SWEEPSET_SIZE(i),bag);
    PetscBagGetData(*bag,(void**)&sweep);
    
    rewind(fp);
    for(i=0; fscanf(fp,"%f",&value) != EOF; i++ )
      sweep->param[i]=(PetscReal)value;
    
    PetscFClose(PETSC_COMM_SELF,fp);
  }
  else
  { 
    i=1; 
    PetscBagCreate(PETSC_COMM_WORLD,SWEEPSET_SIZE(i),bag);
    PetscBagGetData(*bag,(void**)&sweep);
    sweep->param[0]=1;
  }
  PetscBagRegisterInt(*bag,&sweep->size,i,"Size","Number of parameter sweep values");
  PetscBagRegisterRealArray(*bag,&sweep->param,sweep->size,"Parameters","Parameter sweep values");
}

void loadMatricies( const char *optStringArray[], BaseMat baseMatrixArray[], const PetscInt baseMatrixArraySize )
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
  return;
}
