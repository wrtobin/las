#ifndef LAS_NNZ_H_
#define LAS_NNZ_H_
#include <vector>
namespace las
{
  /**
   * Utility struct to store the number of nonzeros per local row in the
   *  diagonal portion of the matrix (dnnz) and the off-diagonal portion
   *  of the matrix. Cast to Sparsity* and pass to PETSc matrix creation
   *  functions, etc.
   *  TODO: build these from SCOREC/core data structures (mesh/field/numbering)
   */
  struct NNZ
  {
    std::vector<int> dnnz;
    std::vector<int> onnz;
  };
}
#endif
