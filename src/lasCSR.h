#ifndef LAS_CSR_H_
#define LAS_CSR_H_
#include "lasSys.h"
#include <vector>
#include <iostream>
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
    /**
     * @param r The number of rows in the sparse matrix
     * @param c The number of cols in the sparse matrix
     * @param nnz The number of nonzeros in the sparse matrix
     * @param rs An array of length r+1 containing the row offsets into cs
     * @param cs An array of length nnz containing the column ids for each nonzero in the matrix
     * @note See the lasCSRBuilder.h file for a basic interface to build these array from the
     *       set of (row,col) values of all nonzeros in a matrix.
     * @note The rs and cs arrays should use 1-indexing for fortran interoperability, if the
     *       first rs offset is zero, it is assumed all values in rs and cs are 0-indexed and
     *       they are converted to use 1-indexing (in a debug build this generates a warning).
     */
    CSR(int r, int c, int nnz, int * rs, int * cs);
    int getNumRows() const { return nr; }
    int getNumCols() const { return nc; }
    int getNumNonzero() const { return nnz; }
    int operator()(int rw, int cl) const
    {
      int result = -1;
      int fst = rws[rw] - 1;
      while((fst < rws[rw+1] - 2) && (cls[fst] - 1 < cl))
        ++fst;
      // the column is correct at offset and the row isn't empty
      if(cls[fst] - 1 == cl && rws[rw] - 1 <= rws[rw+1] - 2)
        result = fst;
      else
        result = -1;
      return result;
    }
    int * getRows() { return &rws[0]; }
    int * getCols() { return &cls[0]; }
    int getMaxEntPerRow();
    int getMaxEntPerCol();
  };
  /*
   * Construct a csr sparse matrix structure from csr arrays
   * \param rws number of rows int the csr matrix
   * \param cls number of columns in the csr matrix
   * \param nnz number of nonzeroes in the csr matrix
   * \param row_arr array where first entry is 1 and all successive entries
   *        hold row_arry[i-1]+number of entries on i-1 row. (one indexed)
   * \param col_arr column index of each matrix entry (should be len nnz)
   */ 
  Sparsity * csrFromArray(int rws, int cls, int nnz, int * row_arr, int * col_arr);
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
