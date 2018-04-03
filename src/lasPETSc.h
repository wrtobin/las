#ifndef LAS_PETSC_H_
#define LAS_PETSC_H_
#include "las.h"
namespace las
{
  class PetscOps;
  void initPETScLAS(int * argc, char ** argv[], MPI_Comm cm);
  void finalizePETScLAS();
  Mat * createPetscMatrix(unsigned lcl, unsigned gbl, unsigned bs, Sparsity * sprs, MPI_Comm cm);
  Vec * createPetscVector(unsigned lcl, unsigned gbl, unsigned bs, MPI_Comm cm);
  Vec * createLHSVector(Mat * m);
  Vec * createRHSVector(Mat * m);
  void destroyPetscMatrix(Mat * m);
  void destroyPetscVector(Vec * v);
  LasOps<PetscOps> * getPetscOps();
  LasSolve * createPetscLUSolve(MPI_Comm cm);
  LasSolve * createPetscQNSolve(void * a);
  LasMultiply * createPetscMultiply();
  class PetscOps : public LasOps<PetscOps>
  {
  public:
    void _zero(las::Mat * m);
    void _zero(las::Vec * v);
    void _assemble(las::Vec * v, int cnt, int * rws, scalar * vls);
    void _assemble(las::Mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls);
    void _set(las::Vec * v, int cnt, int * rws, scalar * vls);
    void _set(las::Mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls);
    scalar _norm(las::Vec * v);
    scalar _dot(las::Vec * v0, las::Vec * v1);
    void _axpy(scalar a, Vec * x, Vec * y);
    void _get(las::Vec * v, scalar *& vls);
    void _restore(las::Vec * v, scalar *& vls);
  };
  mat_builder getPetscMatBuilder();
  vec_builder getPetscVecBuidler();
}
#include "lasPETSc_impl.h"
#endif
