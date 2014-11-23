#ifndef QEPPS_COMMON
#define QEPPS_COMMON

typedef struct
{
  PetscBool Active;
  Mat Matrix;
} BaseMat;

typedef struct
{
    size_t size;
    PetscReal param[1];
} SweepSet;

#define SWEEPSET_SIZE(x) (sizeof(SweepSet) + (sizeof(PetscReal) * ((x) - 1)))

typedef struct
{
    int num;
    Mat matrix[];
} MatrixComponent;

typedef struct
{
    int num;
    PetscComplex param[];
} ParameterSet;

#define MATRIX_COMPONENT_SIZE(x) ( sizeof(MatrixComponent)+sizeof(Mat)*x )
#define PARAMETER_SET_SIZE(x) ( sizeof(ParameterSet)+sizeof(PetscComplex)*x )

#endif
