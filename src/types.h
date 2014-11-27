#ifndef QEPPS_TYPES
#define QEPPS_TYPES

#include <stdbool.h>
#include <complex.h>
#include <string.h>

#define TO_PETSC_COMPLEX(x)  ( creal(x) + PETSC_i*cimag(x) )
#define TO_DOUBLE_COMPLEX(x) ( PetscRealPart(x) + I*PetscImaginaryPart(x) )

typedef struct
{
    int num;
    Mat matrix[];
} MatrixComponent;

#define MATRIX_COMPONENT_SIZE(x) ( sizeof(MatrixComponent)+sizeof(Mat)*x )

#endif
