#include "slepcpep.h"
#include "common.h"

void incorporateMatrixComponent( Mat A, PetscReal scaleValue, BaseMat Ab, BaseMat Ab_p )
{
  if( Ab.Active )
  {
    if( Ab_p.Active ) {
      MatAXPY( A, scaleValue, Ab.Matrix, DIFFERENT_NONZERO_PATTERN );
    } else {
      MatDuplicate( Ab.Matrix, MAT_COPY_VALUES, &A );
      MatScale( A , scaleValue );
    }
  }
}

