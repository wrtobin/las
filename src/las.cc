#include "las.h"
#include "lasDebug.h"
#include <iostream>
#include <lasPETSc.h>
namespace las
{
  void initLAS(int * argc, char *** args, MPI_Comm cm) {
#ifdef HAVE_PETSC
    initPETScLAS(argc, args, cm);
#endif
  }
  void finalizeLAS() {
#ifdef HAVE_PETSC
    finalizePETScLAS();
#endif
  }
  Vec * LasCreateVec::createRHS(Mat*)
  {
    DBG(std::cerr << "[LAS] : createRHS(Mat*) not implemented for configured backend" << std::endl);
    return nullptr;
  }
  Vec * LasCreateVec::createLHS(Mat*)
  {
    DBG(std::cerr << "[LAS] : createLHS(Mat*) not implemented for configured backend" << std::endl);
    return nullptr;
  }
}
