#ifndef LAS_CSR_IMPL_H_
#define LAS_CSR_IMPL_H_
#include "lasDebug.h"
#include "lasInline.h"
#include <iostream>
#include <algorithm>
#include <cassert>
namespace las
{
  LAS_INLINE CSR::CSR(int r, int c,  int nz, int * rs, int * cs)
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
  LAS_INLINE Sparsity * csrFromArray(int rws, int cls, int nnz, int * row_arr, int * col_arr)
  {
    return reinterpret_cast<Sparsity*>(new CSR(rws,cls,nnz,row_arr,col_arr));
  }
  LAS_INLINE int CSR::getMaxEntPerRow() {
    int maxEntPerRow = 0;
    for(std::size_t i=1; i<rws.size();++i) {
     int entPerRow = rws[i]-rws[i-1];
     if(maxEntPerRow<entPerRow) {
       maxEntPerRow = entPerRow;
     }
    }
    return maxEntPerRow;
  }
  LAS_INLINE int CSR::getMaxEntPerCol() {
    std::vector<int> entPerCol(getNumNonzero(),0); 
    for(std::size_t i=0;i<getNumNonzero();++i) {
      assert((cls[i]) < entPerCol.size());
      assert((cls[i]) >= 0);
      ++entPerCol[cls[i]];
      assert(entPerCol[cls[i]] < getNumCols()+1);
    }
   int maxEntPerCol = *std::max_element(entPerCol.begin(), entPerCol.end());
   return maxEntPerCol;
  }
}
#endif
