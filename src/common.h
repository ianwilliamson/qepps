#ifndef QEPPS_COMMON
#define QEPPS_COMMON

typedef struct
{
  PetscBool Active;
  Mat Matrix;
} BaseMat;

// Idea for the below data type and size computation taken from 
// http://stackoverflow.com/questions/7641698/allocating-struct-with-variable-length-array-member

typedef struct
{
    size_t size;
    PetscReal param[1];
} SweepSet;

#define SWEEPSET_SIZE(x) (sizeof(SweepSet) + (sizeof(PetscReal) * ((x) - 1)))

#endif
