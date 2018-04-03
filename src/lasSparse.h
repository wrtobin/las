#ifndef LAS_SPARSE_H_
#define LAS_SPARSE_H_
#include "las.h"
namespace las
{
  Mat * createCSRMatrix(Sparsity* csr);
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
    void _assemble(Vec * v, int cnt, int * rws, scalar * vls);
    void _assemble(Mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls);
    void _set(Vec * v, int cnt, int * rws, scalar * vls);
    void _set(Mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls);
    scalar _norm(Vec * v);
    scalar _dot(Vec * v0, Vec * v1);
    void _axpy(scalar a, Vec * x, Vec * y);
    void _get(Vec * v, scalar *& vls);
    void _restore(Vec * v, scalar *& vls);
  };
}
#include "lasSparse_impl.h"
#endif
