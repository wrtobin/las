#ifndef LAS_CSR_H_
#define LAS_CSR_H_
#include "lasScalar.h"
#include <vector>
namespace las
{
  class Sparsity;
  /**
   * A utility class to describe the sparse structure of a serial
   *  compressed sparse row (CSR) format matrix,
   *  can be used on in conjunction with any linear storage container
   *  to index a sparse matrix.
   */
  class CSR
  {
  private:
    int nr;
    int nc;
    int nnz;
    std::vector<int> rws;
    std::vector<int> cls;
    CSR();
  public:
    CSR(int r, int c, int nnz, int * rs, int * cs);
    int getNumRows() const { return nr; }
    int getNumCols() const { return nc; }
    int getNumNonzero() const { return nnz; }
    int operator()(int rw, int cl) const
    {
      int result = -1;
      int fst = rws[rw] - 1; // 1-indexing    v also here
      while((fst < rws[rw+1] - 2) && (cls[fst] - 1 < cl))
        ++fst;
      if(cls[fst] - 1 == cl) // 1-indexing
        result = fst;
      else
        result = -1;
      return result;
    }
    int * getRows() { return &rws[0]; }
    int * getCols() { return &cls[0]; }
  };
  /**
   * Construct a CSR sparse matrix structure from a full matrix buffer
   * @param mat Pointer to array of rws x cls scalars containing the full matrix
   * @param rws Number of rows in the full matrix
   * @param cls Number of cls in the full matrix
   */
  Sparsity * csrFromFull(scalar * mat, int rws, int cls);
  /**
   *  Produce a full version of the sparse matrix, this is typically
   *   only used for debugging purposes.
   * @param csr A csr object describing the structure of the sprs_mat
   * @param sprs_mat A buffer of length csr->getNumNonzeros() containing the matrix values
   * @param fll_mat A preallocated buffer of length csr->getNumEqs()^2 which will contain
   *                all nonzero and zero values of the matrix;
   */
  void constructFullMatrix(Sparsity * csr,scalar * sprs_mat,scalar * fll_mat);
}
#include "lasCSR_impl.h"
#endif
