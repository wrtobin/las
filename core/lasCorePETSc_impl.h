#ifndef LAS_CORE_PETSC_IMPL_H_
#define LAS_CORE_PETSC_IMPL_H_
#include "lasSparseCore.h"
#include <algorithm>
#include <apfShape.h>
#include <PCU.h>
namespace las
{
  template <>
  Sparsity * createSparsity<petsc>(apf::Numbering * num,
                                   MPI_Comm cm,
                                   int bld_id)
  {
    (void)bld_id;
    return createNNZ(num,cm);
  }
}
#endif
