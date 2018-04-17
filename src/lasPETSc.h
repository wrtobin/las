#ifndef LAS_PETSC_H_
#define LAS_PETSC_H_
#include "las.h"
namespace las
{
  class petsc;
  void initPETScLAS(int * argc, char ** argv[], MPI_Comm cm);
  void finalizePETScLAS();
  template <>
  LasCreateMat * getMatBuilder<petsc>(int id);
  template <>
  LasCreateVec * getVecBuilder<petsc>(int id);
  template <>
  LasOps<petsc> * getLASOps();
  Solve * createPetscLUSolve(MPI_Comm cm);
  Solve * createPetscQNSolve(void * a);
  MatVecMult * createPetscMatVecMult();
  MatMatMult * createPetscMatMatMult();
}
#include "lasPETSc_impl.h"
#endif
