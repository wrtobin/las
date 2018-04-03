#ifndef LAS_SPARSKIT_H_
#define LAS_SPARSKIT_H_
#include "lasSparskitExterns.h"
#include "lasSparse.h"
namespace las
{
  typedef csrOps skOps;
  Mat * createSparskitMatrix(Sparsity * csr);
  void destroySparskitMatrix(Mat * m);
  Vec * createSparskitVector(unsigned n);
  void destroySparskitVector(Vec * v);
  LasOps<skOps> * initSparskitOps();
  /**
   * @note Sparskit is a single-precision solver, so don't use this for anything nonlinear where you're
   *  trying to converge past 10e-8
   */
  LasSolve * createSparskitLUSolve(SparskitBuffers * b);
  LasSolve * createSparskitQuickLUSolve(SparskitBuffers * b);
  LasSolve * createSparskitQuickLUSolve(LasSolve * slvr);
  void printSparskitMat(std::ostream &, Mat * m);
  double getSparskitMatValue(Mat *, int rr, int cc);
  void setSparskitMatValue(Mat *, int rr, int cc, double vl);
  mat_builder getSparskitMatBuilder();
  vec_builder getSparskitVecBuilder();
}
#include "lasSparskit_impl.h"
#endif
