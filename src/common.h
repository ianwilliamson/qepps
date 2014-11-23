#ifndef QEPPS_COMMON
#define QEPPS_COMMON

#define TO_PETSC_COMPLEX(x) ( creal(x)+PETSC_i*cimag(x) )
#define TO_COMPLX_DOUBLE(x) ( PetscReal(x)+I*PetscImag(x) )

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
