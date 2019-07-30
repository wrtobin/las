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
  template <>
  MatDiagonal * getMatDiagonal<sparse>();
  template <>
  MatDiagonalInverse * getMatDiagonalInverse<sparse>();
  template <>
  HadamardProduct * getHadamardProduct<sparse>();
  template <>
  void finalizeMatrix<sparse>(Mat * mat);
  template <>
  void finalizeVector<sparse>(Vec * vec);
  template <>
  void destroySparsity<sparse>(Sparsity * sprs);
}
#include "lasSparse_impl.h"
#endif
