#include "lasCorePETSc.h"
#include <PCU.h>
#include <apf.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <cassert>
#include <iostream>
#include "lasNNZ.h"
namespace las
{
  Sparsity * createPetscSparsity(apf::Numbering * num,
                                 unsigned ndofs,
                                 MPI_Comm cm)
  {
    // number all the nodes
    apf::MeshIterator * it = nullptr;
    apf::MeshEntity * ent = nullptr;
    apf::Mesh * msh = apf::getMesh(num);
    apf::Field * fld = apf::getField(num);
    unsigned blk_sz = apf::countComponents(fld);
    int lmt_dim = -1;
    for (lmt_dim = 3; lmt_dim >= 0; --lmt_dim)
      if (msh->count(lmt_dim) != 0) break;
    apf::Adjacent adj;
    // this may not be correct if we are running in parallel!
    int lcl_blks = ndofs / blk_sz;
    int frst_lcl_blk = 0;
    MPI_Exscan(&lcl_blks, &frst_lcl_blk, 1, MPI_INTEGER, MPI_SUM, cm);
    apf::setNumberingOffset(num, frst_lcl_blk);  // add # to all owned nodes
    apf::synchronize(num);  // get ghost node #s from owners
    NNZ * nnz = new NNZ;
    // build nonzero structure of the matrix
    nnz->dnnz.resize(lcl_blks, 0);
    nnz->onnz.resize(lcl_blks, 0);
    PCU_Comm_Begin();
    it = msh->begin(0);
    while ((ent = msh->iterate(it)))  // only vertices hold nodes
    {
      // FIXME...We assume that the numbering is a naive ordering
      int blk_id = apf::getNumber(num, ent, 0, 0) / blk_sz;
      int adj_cnts[2] = {0, 0};
      // get nodes adjacent to the current node via the
      // element dimension
      apf::getBridgeAdjacent(msh, ent, lmt_dim, 0, adj);
      for (size_t ii = 0; ii < adj.getSize(); ++ii)
        adj_cnts[msh->isOwned(adj[ii])]++;
      int adj_own = adj_cnts[1];
      int adj_gst = adj_cnts[0];
      if (msh->isOwned(ent))
      {
        // assuming 1 node per vtx
        int lcl_blk_id = blk_id - frst_lcl_blk;
        nnz->dnnz[lcl_blk_id] += adj_own + 1;  // adj plus diag blk
        nnz->onnz[lcl_blk_id] += adj_gst;
      }
      else
      {
        int gst_dnnz = 0;
        int ownr = msh->getOwner(ent);
        for (size_t vrt = 0; vrt < adj.getSize(); ++vrt)
        {
          int adj_ownr = msh->getOwner(adj[vrt]);
          // if the adj ghost node is also owned by the same ownr as the primary
          // nd being looped over
          if (adj_ownr == ownr)
          {
            apf::Adjacent edg_frm_frst;
            apf::Adjacent edg_frm_scnd;
            msh->getAdjacent(ent, 1, edg_frm_frst);
            msh->getAdjacent(adj[vrt], 1, edg_frm_scnd);
            for (size_t e1 = 0; e1 < edg_frm_frst.getSize(); ++e1)
              for (size_t e2 = 0; e2 < edg_frm_scnd.getSize(); ++e2)
                if (edg_frm_frst[e1] == edg_frm_scnd[e2])
                  if (msh->isOwned(edg_frm_frst[e1])) gst_dnnz++;
          }
        }
        int adj_lcl = adj_own + adj_gst - gst_dnnz;
        PCU_COMM_PACK(ownr, blk_id);
        PCU_COMM_PACK(ownr, gst_dnnz);
        PCU_COMM_PACK(ownr, adj_lcl);
      }
    }
    PCU_Comm_Send();
    while (PCU_Comm_Receive())
    {
      int blk_id = -1;
      int num_dnnz = 0;
      int num_onnz = 0;
      PCU_COMM_UNPACK(blk_id);
      PCU_COMM_UNPACK(num_dnnz);
      PCU_COMM_UNPACK(num_onnz);
      // this is valid becuase the blk_id is related to an OWNED node
      int lcl_blk_id = blk_id - frst_lcl_blk;
      nnz->dnnz[lcl_blk_id] += num_dnnz;
      nnz->onnz[lcl_blk_id] += num_onnz;
    }
    nnz->rws_size = ndofs;
    nnz->cls_size = ndofs;
    nnz->blk_sz = blk_sz;
    msh->end(it);
    Sparsity * sprs = reinterpret_cast<Sparsity *>(nnz);
    return sprs;
  }
}  // namespace las
