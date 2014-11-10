#ifndef QEPPS_LOAD
#define QEPPS_LOAD

#include "slepcqep.h"
#include "const_qepps.h"

PetscErrorCode loadSweepParameters( PetscInt *nParams, PetscComplex vec_params[] );
PetscErrorCode loadMatricies( const char *optStringArray[], BaseMat baseMatrixArray[], const PetscInt baseMatrixArraySize );

#endif
