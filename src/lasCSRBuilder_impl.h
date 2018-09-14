#ifndef LAS_CSR_BUILDER_IMPL_H_
#define LAS_CSR_BUILDER_IMPL_H_
#include "lasInline.h"
#include "lasCSR.h"
#include <cassert>
#include <iostream>
#include <algorithm>
namespace las
{
  template <typename T1, typename T2>
  bool sort_comps(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
    bool row_less = left.first < right.first;
    bool row_equal = left.first == right.first;
    bool col_less = left.second < right.second;
    return (row_less || (row_equal && col_less));
  }
  template <typename T1, typename T2>
  bool unique_comps(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
    return (left.first == right.first) && (left.second==right.second);
  }
  LAS_INLINE Sparsity * CSRBuilder::finalize()
  {
    
    // we only need to set this size once we know the number of nonzeros!
    // put rows and columns in correct order for adding to the
    // csr structures
    std::sort(coords.begin(), coords.end(), sort_comps<int, int>);
    // get rid of any unique pairs to avoid problems with incorrectly
    auto it = std::unique(coords.begin(), coords.end(), unique_comps<int,int>);
    coords.resize(std::distance(coords.begin(), it));
    assert(nnz >= coords.size());
    if(coords.size() < nnz) {
      std::cerr<<"Warning: ignored "<<nnz-coords.size() << " duplicate entries\n";
    }
    nnz = coords.size();
    cls.resize(nnz);
    for(std::size_t i=0; i<coords.size(); ++i) {
      ++rws[coords[i].first+1];
      cls[i]+=coords[i].second+1;
    }
    for(std::size_t i=1; i<rws.size(); ++i) {
      rws[i] = rws[i]+rws[i-1]-1;
    }
    Sparsity * rslt = reinterpret_cast<Sparsity*>(new CSR(rw_cnt,cl_cnt,nnz,&rws[0],&cls[0]));
    // we may want to pull this out since we don't always want to reuse the same csrbuilder.
    // for now I am leaving this here to avoid any subtle bugs that I might miss.
    reset(); 
    return rslt;
  }
  LAS_INLINE void CSRBuilder::reset()
  {
    coords.clear();
    std::fill(rws.begin(),rws.end(),1);
    nnz = 0;
  }
  LAS_INLINE bool CSRBuilder::need(int rw, int cl)
  {
    bool result = true;
    // negative rw and cl correspond to fixed dofs and won't add a nonzero
    if(rw < 0 || cl < 0)
      result = false;
    if(rw>=rw_cnt || cl>=cl_cnt) 
    {
      std::cerr<<"Warning: inserting a row outside the matrix bounds. Skipping it.\n";
      result=false;
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
    bool result = need(rw,cl);
    if(result)
    {
      coords.push_back(std::pair<int,int>(rw,cl));
      nnz++;
    }
    return result;
  }
}
#endif
