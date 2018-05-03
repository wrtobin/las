#ifndef LAS_SPARSE_DENSE_H_
#define LAS_SPARSE_DENSE_H_
#include "lasSparse.h"
#include "lasDense.h"
namespace las
{
  /**
   * Create a matrix multiplier object which multiplies
   *  a CSR sparse matrix and a dense matrix together and
   *  creates a dense matrix with the result.
   */
  MatMatMult * getSparseMatDenseMatMult();
}
#endif
