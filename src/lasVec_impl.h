#ifndef LAS_VEC_IMPL_H_
#define LAS_VEC_IMPL_H_
#include "lasAlloc.h"
#include "lasInline.h"
#include <lasSys.h>
#include <cassert>
#include <cstring> // memset
#include <iostream>
namespace las
{
  class lasVec
  {
  private:
    scalar * vls;
    int cnt;
    // this is a data location that any negative indices will 
    // sum into
    scalar dummy_data;
  public:
    lasVec(int n)
      : vls(nullptr)
      , cnt(n)
      , dummy_data(0)
    {
      alloc<Malloc>((void**)&vls,sizeof(scalar)*n);
    }
    lasVec(scalar * data, int n)
      : vls(data)
      , cnt(n)
      , dummy_data(0)
    {
    }
    ~lasVec()
    {
      dealloc<Malloc>((void**)vls);
    }
    scalar & operator[](int idx)
    {
      assert(idx < cnt);
      if(idx < 0)
        return dummy_data;
      return vls[idx];
    }
    void setVls(scalar * values) {
      dealloc<Malloc>((void**)vls);
      vls = values;
    }
    int size()
    {
      return cnt;
    }
    // too coupled to the implementation to leave external
    void zero()
    {
      memset(&vls[0],0,sizeof(scalar)*cnt);
      dummy_data = 0;
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
  LAS_INLINE Vec * createVector(scalar * data, unsigned n)
  {
    return reinterpret_cast<Vec*>(new lasVec(data,n));
  }
  LAS_INLINE void destroyVector(Vec * v)
  {
    delete getLASVec(v);
  }
}
#endif
