#ifndef LAS_DEBUG_H_
#define LAS_DEBUG_H_
// move to lasSys.h?
// might not be the correct level of abstraction...
/// las sys is currently for
#ifndef NDEBUG
#define DBG(X) X
#include <cassert>
#else
#define DBG(X)
#endif
#endif
