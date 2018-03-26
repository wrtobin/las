#ifndef LAS_DEBUG_H_
#define LAS_DEBUG_H_
#ifndef NDEBUG
#define DBG(X) X
#include <cassert>
#else
#define DBG(X)
#endif
#endif
