#ifndef LAS_CORE_PETSC_H_
#define LAS_CORE_PETSC_H_
#include <apfNumbering.f>
namespace las
{
  Sparsity * createNNZ(apf::Numbering * num, MPI_Comm cm);
  Mat * createPetscMatrix(apf::Numbering * num, bool owned, MPI_Comm cm);
}
#endif
