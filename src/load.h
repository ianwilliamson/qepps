#ifndef QEPPS_LOAD
#define QEPPS_LOAD

PetscErrorCode loadSweepParameters( PetscInt *nParams, PetscComplex vec_params[] );
PetscErrorCode loadMatricies( const char *optStringArray[], BaseMat baseMatrixArray[], const PetscInt baseMatrixArraySize );

#endif
