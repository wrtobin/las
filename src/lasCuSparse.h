#ifndef LAS_CUSPARSE_H_
#define LAS_CUSPARSE_H_
#include "lasSparse.h"
#include "lasAlloc.h"
namespace las
{
  typedef csrOps cuOps;
  LasOps<cuOps> * getCuSparseOps();
  Mat * createCuMat(CSR * csr);
  void destroyCuMat(Mat * m);
  Vec * createCuVec(unsigned n);
  void destroyCuVec(Vec * v);
  LasSolve * createCuSparseSolve();
  LasMultiply * createCuSparseMultiply();
}
#include "lasCuSparse_impl.h"
#endif
