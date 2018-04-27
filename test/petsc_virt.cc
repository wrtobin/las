#include "petsc_virt.h"
#include <petscmat.h>
class petsc_ops : public ops
{
  void add(las::Mat * m, int rcnt, int * rnum, int ccnt, int * cnum, double * vals)
  {
    MatSetValuesBlocked(*(reinterpret_cast<Mat*>(m)),rcnt,rnum,ccnt,cnum,vals,ADD_VALUES);
  }
};
ops * createPetscOps()
{
  return new petsc_ops;
}
