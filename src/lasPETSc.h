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
  Solve * createPetscQuickLUSolve(Solve * slvr);
  Solve * createPetscQNSolve(void * a);
  MatVecMult * createPetscMatVecMult();
  MatMatMult * createPetscMatMatMult();
  ScalarMatMult * createPetscScalarMatMult();
  template <>
  void finalizeMatrix<petsc>(Mat * mat);
  template <>
  void finalizeVector<petsc>(Vec * vec);
  las::Mat * createPetscMatrix(unsigned l,
                               unsigned bs = 1,
                               Sparsity * sprs = nullptr,
                               MPI_Comm cm = LAS_COMM_WORLD);
  template <>
  void destroySparsity<petsc>(Sparsity * sprs);
}
#include "lasPETSc_impl.h"
#endif
