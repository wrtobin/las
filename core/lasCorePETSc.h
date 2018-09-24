#ifndef LAS_CORE_PETSC_H_
#define LAS_CORE_PETSC_H_
#include "las.h"
#include <apfNumbering.h>
namespace las
{
  //Mat * createPetscMatrix(apf::Numbering * num, bool owned, MPI_Comm cm);
  Sparsity * createPetscSparsity(apf::Numbering * num, unsigned ndofs, MPI_Comm cm=LAS_COMM_WORLD);
}
#endif
