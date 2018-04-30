#ifndef LAS_CSR_BUILDER_IMPL_H_
#define LAS_CSR_BUILDER_IMPL_H_
#include "lasInline.h"
#include "lasCSR.h"
namespace las
{
  LAS_INLINE Sparsity * CSRBuilder::finalize()
  {
    Sparsity * rslt = reinterpret_cast<Sparsity*>(new CSR(rw_cnt,cl_cnt,nnz,&rws[0],&cls[0]));
    reset();
    return rslt;
  }
  LAS_INLINE void CSRBuilder::reset()
  {
    std::fill(rws.begin(),rws.end(),1);
    std::fill(cls.begin(),cls.end(),0);
    nnz = 0;
  }
  LAS_INLINE bool CSRBuilder::need(int rw, int cl)
  {
    bool result = true;
    // negative rw and cl correspond to fixed dofs and won't add a nonzero
    if(rw < 0 || cl < 0)
      result = false;
    else
    {
      for(int ii = rws[rw]; ii < rws[rw+1]; ii++)
      {
        if(cls[ii-1] == (cl+1))
        {
          result = false;
          break;
        }
      }
    }
    return result;
  }
  LAS_INLINE int CSRBuilder::find(int rw, int cl)
  {
    int loc = 0;
    for(loc = rws[rw]; (loc < rws[rw+1]) && (cls[loc-1] < cl+1); ++loc) {}
    return loc;
  }
  LAS_INLINE bool CSRBuilder::add(int rw, int cl)
  {
    bool result = false;
    if((result = need(rw,cl)))
    {
      int ii = find(rw,cl) - 1;
      for(int jj = rw+1; jj < rw_cnt+1; ++jj)
        rws[jj]++;
      for(int jj = rws[rw_cnt] - 2; jj > ii; --jj)
        cls[jj] = cls[jj-1];
      cls[ii] = cl+1;
      nnz++;
    }
    return result;
  }
}
#endif
