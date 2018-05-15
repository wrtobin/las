#ifndef LAS_PETSC_VIRT_H_
#define LAS_PETSC_VIRT_H_
#include <las.h>
#include <petscmat.h>
// basic function call
void add(Mat m, int rcnt, int * rnum, int ccnt, int * cnum, double * vals);

// typedef for c-struct
typedef void (*add_func)(Mat,int,int*,int,int*,double*);
// c-style inheritance through function pointer
//  the actual function is not declared in the header
struct cops
{
  add_func add;
};
// construct a c-style struct api
cops * createPetscCops();

// pure-virtual C++ API
class ops
{
public:
  virtual void add(las::Mat * m, int rcnt, int * rnum, int ccnt, int * cnum, double * vals) = 0;
};
// main is only aware of the super-type
ops * createPetscOps();

#endif
