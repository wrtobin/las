#ifndef LAS_CUSPARSE_H_
#define LAS_CUSPARSE_H_
#include "lasSparse.h"
#include "lasAlloc.h"
namespace las
{
  typedef csrOps cusparse;
  LasOps<cusparse> * getCuSparseOps();
  Mat * createCuMat(CSR * csr);
  void destroyCuMat(Mat * m);
  Vec * createCuVec(unsigned n);
  void destroyCuVec(Vec * v);
  Solve * createCuSparseSolve();
  MatVecMult * createCuMatVecMult();
  MatMatMult * createCuCsrMatCsrMatMult();
  MatMatMult * createCuCsrMatDMatMult();
}
#include "lasCuSparse_impl.h"
#endif
