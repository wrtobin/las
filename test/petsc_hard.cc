#include "petsc_virt.h"
#include "timer.h"
#include <petsc.h>
#include <mpi.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <cassert>
#include <iostream>
#include <string>
int main(int argc, char * argv[])
{
  assert(argc == 3);
  MPI_Init(&argc,&argv);
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  apf::Mesh * msh = apf::loadMdsMesh(argv[1],argv[2]);
  int lmt_dim = msh->getDimension();
  apf::Field * fld = apf::createLagrangeField(msh,"u",apf::VECTOR,1);
  apf::FieldShape * shp = apf::getShape(fld);
  apf::Numbering * num = apf::createNumbering(msh,"u_num",shp,1); // only number the nodes
  apf::MeshIterator * it = nullptr;
  apf::MeshEntity * ent = nullptr;
  int lcl_blks = 0;
  for(int dd = 0; dd < lmt_dim; ++dd)
  {
    if(shp->hasNodesIn(dd))
    {
      it = msh->begin(dd);
      while((ent = msh->iterate(it)))
      {
        if(msh->isOwned(ent))
        {
          int ent_tp = msh->getType(ent);
          int tp_nds = shp->countNodesOn(ent_tp);
          for(int nd = 0; nd < tp_nds; ++nd) // assuming 1 component
            apf::number(num,ent,nd,0,lcl_blks++);
        }
      }
    }
  }
  std::vector<int> dnnz(lcl_blks,0);
  std::vector<int> onnz(lcl_blks,0);
  apf::Adjacent adj;
  it = msh->begin(0);
  while((ent = msh->iterate(it))) // only vertices hold nodes
  {
    if(msh->isOwned(ent))
    {
      //assuming 1 node per vert
      int blk_id = apf::getNumber(num,ent,0,0);
      apf::getBridgeAdjacent(msh,ent,lmt_dim,0,adj);
      for(size_t ii = 0; ii < adj.getSize(); ++ii)
      {
        if(msh->isOwned(adj[ii]))
          dnnz[blk_id]++;
        else
          onnz[blk_id]++;
      }
    }
  }
  int frst_lcl_blk = 0;
  MPI_Exscan(&lcl_blks,&frst_lcl_blk,1,MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD);
  apf::setNumberingOffset(num,frst_lcl_blk);
  apf::synchronize(num); // number ghosts
  int blk_sz = apf::countComponents(fld);
  int lcl_dof_cnt = blk_sz * lcl_blks;
  Mat K;
  MatCreate(PETSC_COMM_WORLD,&K);
  MatSetType(K,MATMPIBAIJ);
  MatSetSizes(K,lcl_dof_cnt,lcl_dof_cnt,PETSC_DETERMINE,PETSC_DETERMINE);
  MatSetBlockSize(K,blk_sz);
  msh->end(it);
  MatMPIBAIJSetPreallocation(K,blk_sz,0,&dnnz[0],0,&onnz[0]);
  MatSetOption(K,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
  it = msh->begin(lmt_dim); // assuming all elements of the same type
  ent = msh->iterate(it);
  apf::MeshElement * msh_lmt = apf::createMeshElement(msh,ent);
  apf::Element * lmt = apf::createElement(fld,msh_lmt);
  int nds_per_lmt = apf::countNodes(lmt);
  int dofs_per_lmt  = (nds_per_lmt * blk_sz) * (nds_per_lmt * blk_sz);
  std::vector<double> ke(dofs_per_lmt,1.0);
  apf::destroyElement(lmt);
  apf::destroyMeshElement(msh_lmt);
  msh->end(it);
  it = msh->begin(lmt_dim);
  apf::NewArray<int> nums(nds_per_lmt);
#if defined(TEST_LASOPS)
  las::lasOps<las::petsc> * las_ops = las::getLasOps<las::petsc>();
  las::Mat * las_K = reinterpret_cast<las::Mat*>(K);
#elif defined(TEST_VIRTUAL)
  ops * petsc_ops = createPetscOps();
#endif
  MatZeroEntries(K); // collective on petsc_comm_world
  unsigned long long span[2] = {0,0};
  span[0] = rdtsc();
#if defined(TEST_RAW)
  while((ent = msh->iterate(it)))
  {
    apf::getElementNumbers(num,ent,nums);
    MatSetValues(K,nds_per_lmt,&nums[0],nds_per_lmt,&nums[0],&ke[0],ADD_VALUES);
  }
  msh->end(it);
#elif defined(TEST_LASOPS)
  while((ent = msh->iterate(it)))
  {
    apf::getElementNumbers(num,ent,nums);
    las_ops->assemble(las_K,nds_per_lmt,&nums[0],nds_per_lmt,&nums[0],&ke[0]);
  }
  msh->end(it);
#elif defined(TEST_VIRTUAL)
  while((ent = msh->iterate(it)))
  {
    apf::getElementNumbers(num,ent,nums);
    petsc_ops->add(las_K,nds_pet_lmt,&nums[0],nds_pet_lmt,&nums[0],&ke[0]);
  }
  msh->end(it);
#endif
  span[1] = rdtsc();
  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);
  MatDestroy(&K);
#if defined(TEST_RAW)
  std::string fnm("raw_results");
#elif defined(TEST_LASOPS)
  std::string fnm("las_results");
#elif defined(TEST_VIRTUAL)
  std::string fnm("virtual_results");
#else
  std::string fnm("oops");
#endif
  int rnk = -1;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rnk);
  int per_rnk_offset = sizeof(unsigned long long) * 2;
  int lcl_offset = rnk * per_rnk_offset;
  MPI_Datatype out_tp;
  MPI_Type_contiguous(2,MPI_UNSIGNED_LONG_LONG,&out_tp);
  MPI_Type_commit(&out_tp);
  MPI_File fout;
  MPI_File_open(PETSC_COMM_WORLD,fnm.c_str(),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fout);
  MPI_File_set_view(fout,lcl_offset,MPI_UNSIGNED_LONG_LONG,out_tp,"native",MPI_INFO_NULL);
  MPI_Status sts;
  MPI_File_write(fout,&span[0],1,out_tp,&sts);
  MPI_File_close(&fout);
  PetscFinalize();
  MPI_Finalize();
}
