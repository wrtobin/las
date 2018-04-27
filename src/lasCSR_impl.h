#ifndef LAS_CSR_IMPL_H_
#define LAS_CSR_IMPL_H_
#include "lasDebug.h"
#include <iostream>
namespace las
{
  inline CSR::CSR(int r, int c,  int nz, int * rs, int * cs)
    : nr(r)
    , nc(c)
    , nnz(nz)
    , rws(rs,rs+nr+1)
    , cls(cs,cs+nz)
  {
    if(rws[0] == 0)
    {
      DBG(std::cerr << "CSR row matrix begins with 0, but the las::CSR type uses 1-indexing, all row and col values will be converted." << std::endl);
      for(int fst = 0; fst < nr+1; ++fst)
        rws[fst]++;
      for(int cl = 0; cl < nz; ++cl)
        cls[cl]++;
    }
  }
  inline Sparsity * csrFromArray(int rws, int cls, int nnz, int * row_arr, int * col_arr)
  {
    return reinterpret_cast<Sparsity*>(new CSR(rws,cls,nnz,row_arr,col_arr));
  }
}
#endif
