#include "lasCore.h"
#include "las.h"
namespace las
{
  int countScopedDofs(apf::Numbering * num, MPI_Comm cm, apf::Sharing * shr = NULL);
  template <class BAK>
  Sparsity * createSparsity(apf::Numbering * num, MPI_Comm cm, int bld_id);
}
