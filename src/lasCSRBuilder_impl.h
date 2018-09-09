#ifndef LAS_CSR_BUILDER_IMPL_H_
#define LAS_CSR_BUILDER_IMPL_H_
#include "lasInline.h"
#include "lasCSR.h"
#include <cassert>
namespace las
{
  LAS_INLINE Sparsity * CSRBuilder::finalize()
  {
    for(std::size_t i=0; i<coords.size(); ++i) {
      ++rws[coords[i].first];
      cls[i] = coords[i].second;
    }
    for(std::size_t i=2; i<rws.size()+1; ++i) {
      rws[i] = rws[i]+rws[i-1]-1;
    }
    Sparsity * rslt = reinterpret_cast<Sparsity*>(new CSR(rw_cnt,cl_cnt,nnz,&rws[0],&cls[0]));
    reset();
    return rslt;
  }
  LAS_INLINE void CSRBuilder::reset()
  {
    coords.clear();
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
      coords.push_back({rw,cl});
      nnz++;
    }
    return result;
  }
}
#endif
