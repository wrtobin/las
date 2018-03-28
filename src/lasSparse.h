#ifndef LAS_SPARSE_H_
#define LAS_SPARSE_H_
#include "las.h"
#include "lasCSR.h"
namespace las
{
  Mat * createCSRMatrix(CSR * csr);
  void destroyCSRMatrix(Mat * m);
  Vec * createVector(unsigned n);
  void destroyVector(Vec * v);
  class csrOps;
  LasOps<csrOps> * initCSROps();
  class csrOps : public LasOps<csrOps>
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
#include "lasSparse_impl.h"
#endif
