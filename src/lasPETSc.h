#ifndef LAS_PETSC_H_
#define LAS_PETSC_H_
#include "las.h"
namespace las
{
  void initPETScLAS(int * argc, char ** argv[], MPI_Comm cm);
  void finalizePETScLAS();
  Mat * createPetscMatrix(int gbl, int lcl, MPI_Comm cm);
  Vec * createPetscVector(int gbl, int lcl, MPI_Comm cm);
  LasOps * getPetscOps();
  LasSolve * createPetscLUSolve();
  LasSolve * createPetscQNSolve(void * a);
  LasMultiply * createPetscMultiply();
}
#endif
