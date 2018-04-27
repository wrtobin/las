#ifndef LAS_PETSC_VIRT_H_
#define LAS_PETSC_VIRT_H_
#include <las.h>
class ops
{
  virtual void add(las::Mat * m, int rcnt, int * rnum, int ccnt, int * cnum, double * vals) = 0;
};
ops * createPetscOps();
typedef void (*add_func)(Mat,int,int*,int,int*,double*);
struct cops
{
  add_func add;
};
cops * createPetscCops();
void add(Mat m, int rcnt, int * rnum, int ccnt, int * cnum, double * vals);
#endif
