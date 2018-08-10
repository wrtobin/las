#ifndef LAS_CORE_H_
#define LAS_CORE_H_
#include "las.h"
#include <apf.h>
#include <apfAggregateNumbering.h>
namespace las
{
  template <class B>
  Mat * createMatrix(apf::AggNumbering * num,
                     MPI_Comm cm = LAS_COMM_WORLD,
                     int bld_id = 0);
  template <class B>
  Sparsity * createSparsity(apf::AggNumbering * num,
                            MPI_Comm cm = LAS_COMM_WORLD,
                            int bld_id = 0);
  void destroyMatrix(Mat * m);
}
#endif
