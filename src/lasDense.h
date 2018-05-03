#ifndef LAS_DENSE_H_
#define LAS_DENSE_H_
#include "las.h"
namespace las
{
  class dense;
  template <>
  LasCreateMat * getMatBuilder<dense>(int);
  template <>
  LasCreateVec * getVecBuilder<dense>(int);
  template <>
  LasOps<dense> * getLASOps();
  Sparsity * createDensity(int rws, int cls, scalar * vls = nullptr);
  MatVecMult * getDenseMatVecMult();
  MatMatMult * getDenseMatMatMult();
}
#include "lasDense_impl.h"
#endif
