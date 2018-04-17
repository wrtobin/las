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
}
#include "lasSparse_impl.h"
#endif
