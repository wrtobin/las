#ifndef LAS_SPARSKIT_H_
#define LAS_SPARSKIT_H_
#include "lasSparskitExterns.h"
#include "lasSparse.h"
#include <ostream>
namespace las
{
  typedef sparse sparskit;
  template <>
  LasCreateMat * getMatBuilder<sparskit>(int id);
  template <>
  LasCreateVec * getVecBuilder<sparskit>(int id);
  template <>
  LasOps<sparskit> * getLASOps<sparskit>();
  /**
   * @note Sparskit is a single-precision solver, so don't use this for anything nonlinear where you're
   *  trying to converge past 10e-8
   */
  Solve * createSparskitLUSolve(SparskitBuffers * b);
  Solve * createSparskitQuickLUSolve(SparskitBuffers * b);
  Solve * createSparskitQuickLUSolve(Solve * slvr);
  void printSparskitMat(std::ostream &, Mat * m);
  double getSparskitMatValue(Mat *, int rr, int cc);
  void setSparskitMatValue(Mat *, int rr, int cc, double vl);
}
#include "lasSparskit_impl.h"
#endif
