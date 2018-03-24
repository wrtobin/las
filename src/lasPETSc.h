#ifndef LAS_PETSC_H_
#define LAS_PETSC_H_
#include "las.h"
namespace las
{
  Mat * createPetscMatrix(int gbl, int lcl);
  Vec * createPetscVector(int gbl, int lcl);
  LasOps * getPetscOps();
  LasSolve * createPetscLUSolve();
  LasSolve * createPetscQNSolve(void * a);
  LasMultiply * createPetscMultiply();
}
#endif
