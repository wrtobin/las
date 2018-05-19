#include "petsc_virt.h"
#include "timer.h"
#include <lasConfig.h>
#include <petsc.h>
#include <mpi.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <gmi_sim.h>
#include <PCU.h>
#include <cassert>
#include <iostream>
#include <string>
int main(int argc, char * argv[])
{
  assert(argc == 3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_sim_start();
  gmi_register_sim();
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  //PetscPopSignalHandler();
  int rnk = -1;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rnk);
  apf::Mesh * msh = apf::loadMdsMesh(argv[1],argv[2]);
  int lmt_dim = msh->getDimension();
  apf::Field * fld = apf::createLagrangeField(msh,"u",apf::VECTOR,1);
  apf::zeroField(fld);
  apf::FieldShape * shp = apf::getShape(fld);
  apf::Numbering * num = apf::createNumbering(msh,"u_num",shp,1); // only number the nodes
  // number all the nodes
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
          for(int nd = 0; nd < tp_nds; ++nd) // assuming 1 component / blocks
            apf::number(num,ent,nd,0,lcl_blks++);
        }
      }
      msh->end(it);
    }
  }
  int frst_lcl_blk = 0;
  MPI_Exscan(&lcl_blks,&frst_lcl_blk,1,MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD);
  apf::setNumberingOffset(num,frst_lcl_blk);
  apf::synchronize(num); // number ghosts
  /*
  apf::writeVtkFiles("cube_num",msh);
  int own_vrts = apf::countOwned(msh,0);
  std::cout << "[" << rnk << "] : owns " << own_vrts << " vtx" << std::endl;
  it = msh->begin(0);
  while((ent = msh->iterate(it))) // only vertices hold nodes
    if(!msh->isOwned(ent))
    {
      apf::Copies gsts;
      msh->getRemotes(ent,gsts);
      int gst_cnt = gsts.size();
      std::cout << "[" << rnk << "] : gst vtx with " << gst_cnt << " gsts" << std::endl;
    }
  msh->end(it);
  // report ownership infor for debugging
  it = msh->begin(0);
  while((ent = msh->iterate(it))) // only vertices hold nodes
    if(msh->isOwned(ent))
    {
      int ent_tp = msh->getType(ent);
      int tp_nds = shp->countNodesOn(ent_tp);
      int blk_id = -1;
      for(int nd = 0; nd < tp_nds; ++nd) // assuming 1 component / blocks
      {
        apf::getNumber(num,ent,nd,blk_id);
        std::cout << "[" << rnk << "] : owns vtx with node " << blk_id << std::endl;
      }
    }
  msh->end(it);
  */
  // build nonzero structure of the matrix
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
      int adj_id = -1;
      apf::getBridgeAdjacent(msh,ent,lmt_dim,0,adj);
      for(size_t ii = 0; ii < adj.getSize(); ++ii)
      {
        adj_id = apf::getNumber(num,adj[ii],0,0);
        //std::cout << "[" << rnk << "] : owns " << blk_id;
        if(msh->isOwned(adj[ii]))
        {
          dnnz[blk_id]++;
          //std::cout << " adj to ownd " << adj_id << std::endl;
        }
        else
        {
          onnz[blk_id]++;
          //std::cout << " adj to ghst " << adj_id << std::endl;
        }
      }
    }
  }
  msh->end(it);
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
#if defined(TEST_LASOPS) || defined(TEST_VIRTUAL)
  las::Mat * las_K = reinterpret_cast<las::Mat*>(&K);
#endif
#if defined(TEST_LASOPS)
  las::LasOps<las::petsc> * las_ops = las::getLASOps<las::petsc>();
#elif defined(TEST_CVIRT)
  cops * petsc_cops = createPetscCops();
#elif defined(TEST_VIRTUAL)
  ops * petsc_ops = createPetscOps();
#endif
  MatZeroEntries(K); // collective on petsc_comm_world
#ifdef TEST_SINGLE
  //unsigned long long inst[2] = {0,0};
  double t3 = 0.0;
  double t4 = 0.0;
#else
  //unsigned long long span[2] = {0,0};
    double t1 = PCU_Time();
#endif
  //span[0] = rdtsc();
  while((ent = msh->iterate(it)))
  {
    apf::getElementNumbers(num,ent,nums);
#ifdef TEST_SINGLE
    //inst[0] = rdtsc();
    t3 = PCU_Time();
#endif
#if defined(TEST_RAW)
    MatSetValues(K,nds_per_lmt,&nums[0],nds_per_lmt,&nums[0],&ke[0],ADD_VALUES);
#elif defined(TEST_LASOPS)
    las_ops->assemble(las_K,nds_per_lmt,&nums[0],nds_per_lmt,&nums[0],&ke[0]);
#elif defined(TEST_CALL)
    add(K,nds_per_lmt,&nums[0],nds_per_lmt,&nums[0],&ke[0]);
#elif defined(TEST_CVIRT)
    (*petsc_cops->add)(K,nds_per_lmt,&nums[0],nds_per_lmt,&nums[0],&ke[0]);
#elif defined(TEST_VIRTUAL)
    petsc_ops->add(las_K,nds_per_lmt,&nums[0],nds_per_lmt,&nums[0],&ke[0]);
#endif
#ifdef TEST_SINGLE
    //inst[1] = rdtsc();
    t4 = PCU_Time();
    break;
#endif
  }
  msh->end(it);
#ifndef TEST_SINGLE
  //span[1] = rdtsc();
  double t2 = PCU_Time();
#endif
  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);
  //MatDestroy(&K);
  /*
#if defined(TEST_RAW)
  std::string fnm("raw_results");
#elif defined(TEST_LASOPS)
  std::string fnm("las_results");
#elif defined(TEST_CALL)
  std::string fnm("call_results");
#elif defined(TEST_CVIRT)
  std::string fnm("cvirt_results");
#elif defined(TEST_VIRTUAL)
  std::string fnm("virtual_results");
#else
  std::string fnm("oops");
#endif
  int per_rnk_offset = sizeof(unsigned long long) * 2;
  int lcl_offset = rnk * per_rnk_offset;
  MPI_Datatype out_tp;
  MPI_Type_contiguous(2,MPI_UNSIGNED_LONG_LONG,&out_tp);
  MPI_Type_commit(&out_tp);
  MPI_File fout;
  MPI_File_open(PETSC_COMM_WORLD,fnm.c_str(),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fout);
  MPI_File_set_view(fout,lcl_offset,MPI_UNSIGNED_LONG_LONG,out_tp,"native",MPI_INFO_NULL);
  MPI_Status sts;
#ifdef TEST_SINGLE
    MPI_File_write(fout,&inst[0],1,out_tp,&sts);
#else
    MPI_File_write(fout,&span[0],1,out_tp,&sts);
#endif
  MPI_File_close(&fout);
  */
#ifdef TEST_SINGLE
  //std::cout << CLOCKS_PER_SEC << ", " << inst[0] << ", " << inst[1] << ", " << (double)(inst[1] - inst[0]) / (double)CLOCKS_PER_SEC << std::endl;
  std::cout << "Single API call took " << (t4-t3) << " seconds." << std::endl;
#else
  //std::cout << CLOCKS_PER_SEC << ", " << span[0] << ", " << span[1] << ", " << (double)(span[1] - span[0]) / (double)CLOCKS_PER_SEC << std::endl;
  std::cout << "FEM assembly took " << (t2-t1) << " seconds." << std::endl;
#endif
  PetscFinalize();
  gmi_sim_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}
