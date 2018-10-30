#ifndef LAS_SPARSE_H_
#define LAS_SPARSE_H_
#include "las.h"
namespace las
{
  class sparse;
  template <>
  LasCreateMat * getMatBuilder<sparse>(int);
  template <>
  LasCreateVec * getVecBuilder<sparse>(int);
  template <>
  LasOps<sparse> * getLASOps();
  MatVecMult * getSparseMatVecMult();
  MatMatMult * getSparseMatMatMult();
  ScalarMatMult * getSparseScalarMatMult();
  ScalarMatScalarMatAdd * getSparseScalarMatScalarMatAdd();
  template <>
  void finalizeMatrix<sparse>(Mat * mat);
  template <>
  void finalizeVector<sparse>(Vec * vec);
  template <>
  void destroySparsity<sparse>(Sparsity * sprs);
}
#include "lasSparse_impl.h"
#endif
