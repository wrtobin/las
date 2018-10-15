#include <las.h>
#include <lasPETSc.h>
#include <lasSparskit.h>
#include <mpi.h>
#include <cassert>
#include <iostream>
#include <utility>
#include <vector>
int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  las::initLAS(&argc, &argv);
  auto ops = las::getLASOps<las::BACKEND>();
  auto vb = las::getVecBuilder<las::BACKEND>(0);
  las::Vec * v = vb->create(3, 1, LAS_COMM_WORLD);
  las::Vec * u = vb->create(3, 1, LAS_COMM_WORLD);
  ops->zero(v);
  ops->zero(u);
  double vals[] = {1.0, 2.0, 3.0};
  int rws[] = {0, 1, 2};
  ops->assemble(v, 3, &rws[0], &vals[0]);
  auto mb = las::getMatBuilder<las::BACKEND>(0);
#if defined USING_SPARSKIT
  las::CSRBuilder csrBldr(3, 3);
  for (int i = 0; i < 3; ++i)
  {
    csrBldr.add(rws[i], rws[i]);
  }
  csrBldr.add(rws[0], rws[1]);
  las::Sparsity * sprs = csrBldr.finalize();
#elif defined USING_PETSC
  las::Sparsity * sprs = nullptr;
#endif
  las::Mat * m = mb->create(3, 1, sprs, MPI_COMM_SELF);
  ops->assemble(m, 1, &rws[0], 1, &rws[0], &vals[0]);
  ops->assemble(m, 1, &rws[1], 1, &rws[1], &vals[1]);
  ops->assemble(m, 1, &rws[2], 1, &rws[2], &vals[2]);
  ops->assemble(m, 1, &rws[0], 1, &rws[1], &vals[2]);
  las::finalizeMatrix<las::BACKEND>(m);
  las::finalizeVector<las::BACKEND>(v);
#if defined USING_SPARSKIT
  las::printSparskitMat(std::cout, m, las::PrintType::full);
  int nnz = ((las::CSR *)sprs)->getNumNonzero();
  int nrows = ((las::CSR *)sprs)->getNumRows();
  las::SparskitBuffers bfrs = las::SparskitBuffers(nrows, nnz * 100);
  las::Solve * slv = createSparskitLUSolve(&bfrs, 1E-6);
#elif defined USING_PETSC
  las::Solve * slv = las::createPetscLUSolve(LAS_COMM_WORLD);
#endif
  slv->solve(m, u, v);
  double * rslt = nullptr;
  ops->get(u, 3, &rws[0], &rslt);
  std::cout << "results" << std::endl;
  for (int i = 0; i < 3; ++i)
  {
    std::cout << rslt[i] << " ";
  }
  std::cout << std::endl;
  delete[] rslt;
  mb->destroy(m);
  vb->destroy(u);
  vb->destroy(v);
  delete slv;
  delete mb;
  delete vb;
  delete ops;
  las::destroySparsity<las::BACKEND>(sprs);
//#if BACKEND == sparskit
//  las::destroySparsity<las::CSR *>(sprs);
//#endif
  las::finalizeLAS();
  MPI_Finalize();
}
