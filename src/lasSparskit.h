#ifndef LAS_SPARSKIT_H_
#define LAS_SPARSKIT_H_
#include "lasCSR.h"
#include "lasSparskitExterns.h"
#include <las.h>
namespace las
{
  class skOps;
  Mat * createSparskitMatrix(CSR * csr);
  void deleteSparskitMatrix(Mat * m);
  Vec * createSparskitVector(int n);
  void deleteSparskitVector(Vec * v);
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
  class skOps : public LasOps<skOps>
  {
  public:
    void _zero(Mat * m);
    void _zero(Vec * v);
    void _assemble(Vec * v, int cnt, int * rws, double * vls);
    void _assemble(Mat * m, int cntr, int * rws, int cntc, int * cls, double * vls);
    void _set(Vec * v, int cnt, int * rws, double * vls);
    void _set(Mat * m, int cntr, int * rws, int cntc, int * cls, double * vls);
    double _norm(Vec * v);
    double _dot(Vec * v0, Vec * v1);
    void _axpy(double a, Vec * x, Vec * y);
    void _get(Vec * v, double *& vls);
    void _restore(Vec * v, double *& vls);
  };
}
#include "lasSparskit_impl.h"
#endif
