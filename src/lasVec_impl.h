#ifndef LAS_VEC_IMPL_H_
#define LAS_VEC_IMPL_H_
#include "lasAlloc.h"
#include "lasInline.h"
#include <cassert>
namespace las
{
  class lasVec
  {
  private:
    scalar * vls;
    int cnt;
  public:
    lasVec(int n)
      : vls(nullptr)
      , cnt(n)
    {
      alloc<Malloc>((void**)&vls,sizeof(scalar)*(n+1));
    }
    ~lasVec()
    {
      dealloc<Malloc>((void**)vls);
    }
    scalar & operator[](int idx)
    {
      assert(idx < cnt);
      if(idx < 0)
        idx = cnt;
      return vls[idx];
    }
    int size()
    {
      return cnt;
    }
    // too coupled to the implementation to leave external
    void zero()
    {
      memset(&vls[0],0,sizeof(scalar)*(cnt+1));
    }
  };
  LAS_INLINE lasVec * getLASVec(Vec * v)
  {
    return reinterpret_cast<lasVec*>(v);
  }
  LAS_INLINE Vec * createVector(unsigned n)
  {
    return reinterpret_cast<Vec*>(new lasVec(n));
  }
  LAS_INLINE void destroyVector(Vec * v)
  {
    delete getLASVec(v);
  }
}
#endif
