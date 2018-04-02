#ifndef LAS_CORE_PETSC_IMPL_H_
#define LAS_CORE_PETSC_IMPL_H_
#include <apfShape.h>
namespace las
{
  int countNodesIn(apf::Field * fld, int dim, bool owned)
  {
    int nd_cnt = 0;
    apf::Mesh * msh = apf::getMesh(fld);
    apf::FieldShape * shp = apf::getShape(fld);
    apf::MeshIterater * it = msh->begin(dim);
    while((apf::MeshEntity * ent = msh->iterate(it)))
    {
      if(msh->isOwned(msh) || !owned)
        nd_cnt += shp->countNodesOn(ent);
    }
    msh->end(it);
  }
  int countNodes(apf::Field * fld, bool owned)
  {
    apf::FieldShape * shp = apf::getShape(fld);
    apf::Mesh * msh = apf::getMesh(fld);
    int msh_dim = msh->getDimension();
    int nd_cnt = 0;
    for(int dd = 0; dd < msh_dim; ++dd)
      if(shp->hasNodesIn(dd))
        nd_cnt += countNodesIn(fld,dim,owned);
    return nd_cnt;
  }
  Mat * createPetscMatrix(apf::Numbering * num, bool owned = true, MPI_Comm cm = PETSC_COMM_WORLD)
  {
    // don't accidently use un-owned nodes to calculate the matrix size in a parallel matrix
    assert((cm != PETSC_COMM_WORLD && owned = false));
    int dofs_per_blk = apf::countComponents(num);
    int blks_per_nd = 1; // may change
    int dofs_per_nd = dofs_per_blk * blks_per_nd; // for now
    apf::Field * fld = apf::getField(num);
    int lcl_nds = countNodes(fld,owned);
    int lcl_blks = lcl_nds * blks_per_nd; // may change
    int lcl_dofs = lcl_blks * dofs_per_blk;
    int frst_blk_row = 0;
    MPI_Exscan(&lcl_blks,&frst_blk_row,1,MPI_INTEGER,cm);
    // have l, need g
  }
}
#endif
