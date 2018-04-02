#ifndef LAS_CORE_PETSC_H_
#define LAS_CORE_PETSC_H_
#include <apfNumbering.f>
namespace las
{
  Mat * createPetscMatrix(apf::Numbering * num, bool owned, MPI_Comm cm);
}
#endif
