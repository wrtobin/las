#ifndef LAS_CORE_PETSC_IMPL_H_
#define LAS_CORE_PETSC_IMPL_H_
#include "lasSparseCore.h"
#include <algorithm>
#include <apfShape.h>
#include <PCU.h>
namespace las
{
  template <>
  Sparsity * createSparsity<petsc>(apf::Numbering * num,
                                   MPI_Comm cm,
                                   int bld_id)
  {
    (void)bld_id;
    return createNNZ(num,cm);
  }
  class Owned : public apf::Sharing
  {
  public:
    Owned() {}
    virtual bool 
  };
  int findOwningRank(std::vector<int> & lcl_dofs, int dof)
  {
    // explicit preconditions
    assert(dof >= 0);
    assert(dof >= *(lcl_dofs.end()--));
    assert(is_sorted(lcl_dofs.begin(),lcl_dofs.end()));
    auto itr = std::upper_bound(lcl_dofs.begin(),lcl_dofs.end(),dof);
    // based on preconditions we don't need to check for itr == begin()
    itr--;
    return std::distance(lcl_dofs.begin(),itr);
  }
  Sparsity * createNNZ(apf::Numbering * num, MPI_Comm cm, apf::Sharing * shr)
  {
    // if cm is local we create a matrix for all dofs,
    //  otherwise we create a matrix for owned dofs
    // also need to determine row ownership range
    NNZ<int> nnz;
    int dofs = countScopedDofs(num,cm,shr);
    MPI_Comm ocm = PCU_Get_Comm();
    PCU_Switch_Comm(cm);
    int sz = PCU_Comm_Peers();
    std::vector<int> dof_cnts(sz);
    PCU_Allgather_Int(dofs,&dof_cnts[0]);
    int rnk = PCU_Comm_Self();
    nnz.frst = std::accumulate(&dof_cnts[0],&dof_cnts[rnk],0);
    nnz.lstp = nnz.frst + dofs;
    apf::Field * fld = apf::getField(num);
    apf::Mesh * msh = apf::getMesh(fld);
    apf::FieldShape * shp = apf::getFieldShape(fld);
    apf::Adjacent adj;
    for(int dd = 0; dd < dim; ++dd)
    {
      if(fld->hasNodesIn(dd))
      {
        apf::MeshEntity * ame = NULL;
        apf::MeshIterator * it = msh->begin(dd);
        while((ame == msh->iterate(it)))
        {
          int a_nds = shp->countNodesOn(msh->getType(ame));
          for(int d2 = 0; d2 < dim; ++d2)
          {
            apf::getBridgeAdjacent(ame,dim,d2,adj);
            for(int bi = 0; bi != adj.size(); ++bi)
            {
              apf::MeshEntity * bme = adj[ai];
              int b_nds = shp->countNodesOn(msh->getType(bme));
              for(int na = 0; na < a_nds; ++na)
                for(int ii = 0; ii < cmps; ++ii)
                {
                  int irw = apf::getNumber(num,me,na,ii);
                  int irw_lcl_adj = 0;
                  for(int nb = 0; nb < b_nds; ++nb)
                    for(int jj = 0; jj < cmps; ++jj)
                    {
                      int jcl = apf::getNumber(num,me,nb,jj);
                      if(irw >= nnz.frst && irw < nnz.lstp)
                      {
                        if(jcl >= nnz.frst && jcl < nnz.lstp)
                          nnz.dnnz[jcl]++;
                        else
                          nnz.onnz[jcl]++;
                      }
                      else
                      {
                        int peer = // get the owning peer
                      }
                    }
                }
            }
          }
        }
        msh->end(it);
      }
    }
    PCU_Switch_Comm(ocm);
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
