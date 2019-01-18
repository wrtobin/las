#ifndef LAS_SPARSKIT_H_
#define LAS_SPARSKIT_H_
#include "lasSparskitExterns.h"
#include "lasSparse.h"
#include <ostream>
namespace las
{
  enum class PrintType { full, mmarket };
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
  Solve * createSparskitLUSolve(SparskitBuffers * b, double eps = 0.0);
  Solve * createSpaskitLUSolve(Solve * slvr, double eps = 0.0);
  Solve * createSparskitQuickLUSolve(SparskitBuffers * b, double eps = 0.0);
  Solve * createSparskitQuickLUSolve(Solve * slvr, double eps = 0.0);
  template <>
  MatVecMult * getMatVecMult<sparse>();
  template <>
  MatMatMult * getMatMatMult<sparse>();
  template<>
  ScalarMatMult * getScalarMatMult<sparse>();
  template <>
  MatMatAdd * getMatMatAdd<sparse>();
  template <>
  VecVecAdd * getVecVecAdd<sparse>();
  void printSparskitMat(std::ostream &,
                        Mat * m,
                        PrintType pt = PrintType::full,
                        bool symmetric = false);
  Mat * readSparskitMat(std::istream & istream, PrintType pt = PrintType::full);
  bool sparskitMatClose(Mat * m1, Mat * m2, double rtol=1E-15, double atol=1E-15);
}  // namespace las
#include "lasSparskit_impl.h"
#endif
