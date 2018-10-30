#ifndef LAS_CSRCORE_H_
#define LAS_CSRCORE_H_
#include "lasCSR.h"
#include <apfNumbering.h>
namespace las
{
  /**
   * Construct a CSR sparse matrix structure to manage a sparse matrix based
   *  on the dof numbering of an apf field. This only operates on the local
   *  values of the numbering, no mesh partitioning is considered.
   */
  Sparsity * createCSR(apf::Numbering * num, int ndofs);
  Sparsity * createIdentityCSR(int ndofs);
}
#endif
