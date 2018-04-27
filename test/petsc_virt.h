#ifndef LAS_PETSC_VIRT_H_
#define LAS_PETSC_VIRT_H_
#include <las.h>
class ops
{
  virtual void add(las::Mat * m, int rcnt, int * rnum, int ccnt, int * cnum, double * vals) = 0;
};
ops * createPetscOps();
#endif
