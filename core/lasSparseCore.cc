#include "lasSparseCore.h"
namespace las
{
  int countScopedDofs(apf::Numbering * num, MPI_Comm cm, apf::Sharing * shr)
  {
    apf::Field * fld = apf::getField(num);
    apf::Mesh * msh = apf::getMesh(fld);
    int dim = apf::getDimension(msh);
    bool is_par = cm != MPI_COMM_SELF;
    apf::Sharing * s = is_par ? shr : new apf::AllSharing(msh);
    int nds = 0;
    for(int dd = 0; dd < dim; ++dd)
      nds += apf::countOwnedNodes(msh,dd,apf::getShape(fld),s);
    int cmps = apf::countComponents(fld);
    return nds * cmps;
  }
  template <>
  Sparsity * createSparsity<dense>(apf::Numbering * num, MPI_Comm cm, int bld_id)
  {
    (void)bld_id;
    assert(cm == MPI_COMM_SELF && "createSparsity<dense> only implemented for local matrices!");
    int dofs = countScopedDofs(num,cm);
    return createDensity(dofs,dofs);
  }
  template <>
  Sparsity * createSparsity<sparse>(apf::Numbering * num, MPI_Comm cm, int bld_id)
  {
    (void)bld_id;
    assert(cm == MPI_COMM_SELF && "createSparsity<sparse> only implemented for local matrices!");
    int dofs = countScopedDofs(num,cm);
    return createCSR(num,dofs);
  }
}
