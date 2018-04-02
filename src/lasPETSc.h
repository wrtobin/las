#ifndef LAS_PETSC_H_
#define LAS_PETSC_H_
#include "las.h"
namespace las
{
  class PetscOps;
  void initPETScLAS(int * argc, char ** argv[], MPI_Comm cm);
  void finalizePETScLAS();
  Mat * createPetscMatrix(int gbl, int lcl, int bs, double * dnnz, double * onnz, MPI_Comm cm);
  Vec * createPetscVector(int gbl, int lcl, MPI_Comm cm);
  Vec * createLHSVector(Mat * m);
  Vec * createRHSVector(Mat * m);
  LasOps<PetscOps> * getPetscOps();
  LasSolve * createPetscLUSolve(MPI_Comm cm);
  LasSolve * createPetscQNSolve(void * a);
  LasMultiply * createPetscMultiply();
  class PetscOps : public LasOps<PetscOps>
  {
    void _zero(las::Mat * m);
    void _zero(las::Vec * v);
    void _assemble(las::Vec * v, int cnt, int * rws, double * vls);
    void _assemble(las::Mat * m, int cntr, int * rws, int cntc, int * cls, double * vls);
    void _set(las::Vec * v, int cnt, int * rws, double * vls);
    void _set(las::Mat * m, int cntr, int * rws, int cntc, int * cls, double * vls);
    double _norm(las::Vec * v);
    double _dot(las::Vec * v0, las::Vec * v1);
    void _axpy(double a, Vec * x, Vec * y);
    void _get(las::Vec * v, double *& vls);
    void _restore(las::Vec * v, double *& vls);
  };
}
#include "lasPETSc_impl.h"
#endif
