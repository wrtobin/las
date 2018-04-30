#ifndef LAS_CSR_BUILDER_H_
#define LAS_CSR_BUILDER_H_
#include <vector>
namespace las
{
  class CSR;
  class Sparsity;
  class CSRBuilder
  {
  private:
    std::vector<int> rws;
    std::vector<int> cls;
    int nnz;
    int rw_cnt;
    int cl_cnt;
  public:
    CSRBuilder(int num_rows, int num_cols)
      : rws(num_rows+1,1)
      , cls(num_rows*num_cols,0)
      , nnz(0)
      , rw_cnt(num_rows)
      , cl_cnt(num_cols)
    { }
    bool add(int rw, int cl);
    Sparsity * finalize();
    /**
     * Reset the internal data structure to reuse the
     *  builder for another CSR structure, does not
     *  allow resizing the CSR so only use this when
     *  the dimensions of the matrix are not changing.
     */
    void reset();
  protected:
    bool need(int rw, int cl);
    int find(int rw, int cl);
  };
}
#include "lasCSRBuilder_impl.h"
#endif
