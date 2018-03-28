#ifndef LAS_CSRBUILDER_H_
#define LAS_CSRBUILDER_H_
namespace las
{
  bool needNonzero(int * rws, int rw, int * cls, int cl);
  int findLocation(int * rws, int rw, int * cls, int cl, int neq);
  bool addNonzero(int * rws, int rw, int * cls, int cl, int neq);
}
#endif
