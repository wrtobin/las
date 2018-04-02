#ifndef LAS_ALLOC_H_
#define LAS_ALLOC_H_
#include <cstdlib>
namespace las
{
  template <typename T>
  void alloc(void **, size_t);
  template <typename T>
  void dealloc(void **);
  class Malloc;
  template <>
  inline void alloc<Malloc>(void ** dat, size_t sz)
  {
    (*dat) = malloc(sz);
  }
  template <>
  inline void dealloc<Malloc>(void ** dat)
  {
    free(dat);
  }
}
#endif
