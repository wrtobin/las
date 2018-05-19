#include "petsc_virt.h"
#include <petscmat.h>
// function call implementation
void add(Mat m, int rcnt, int * rnum, int ccnt, int * cnum, double * vals)
{
  MatSetValues(m,rcnt,rnum,ccnt,cnum,vals,ADD_VALUES);
}

// unexposed function used for c-style interface
void _add(Mat m, int rcnt, int * rnum, int ccnt, int * cnum, double * vals)
{
  MatSetValues(m,rcnt,rnum,ccnt,cnum,vals,ADD_VALUES);
}

// create a C-style API struct and bind the function pointer
cops * createPetscCops()
{
  cops * petsc_cops = new cops;
  petsc_cops->add = _add;
  return petsc_cops;
}
// inherit the C++ virtual class
class petsc_ops : public ops
{
  void add(las::Mat * m, int rcnt, int * rnum, int ccnt, int * cnum, double * vals)
  {
    MatSetValues(*(reinterpret_cast<Mat*>(m)),rcnt,rnum,ccnt,cnum,vals,ADD_VALUES);
  }
};
// retrieve the C++ API as the super-type
ops * createPetscOps()
{
  return new petsc_ops;
}
