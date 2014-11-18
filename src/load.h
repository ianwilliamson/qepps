#ifndef QEPPS_LOAD
#define QEPPS_LOAD

void loadSweepParameters( PetscBag *bag );
void loadMatricies( const char *optStringArray[], BaseMat baseMatrixArray[], const PetscInt baseMatrixArraySize );

#endif
