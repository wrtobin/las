#include "lasSparse.h"
#include <cassert>
#include <sstream>
int main(int,char*[])
{
  double eye[] = { 1.0 , 0.0 , 0.0 ,
                   0.0 , 1.0 , 0.0 ,
                   0.0 , 0.0 , 1.0 };
  double col[] = { 0.0 , 3.14, 0.0 ,
                   0.0 , 3.14, 0.0 ,
                   0.0 , 3.14, 0.0 };
  las::Sparsity * eye_csr = las::csrFromFull(&eye[0],3,3);
  las::Sparsity * col_csr = las::csrFromFull(&col[0],3,3);
  assert(eye_csr && col_csr);
  las::LasCreateMat * mat_fct = las::getMatBuilder<las::sparse>(0);
  las::Mat * eye_mat = mat_fct->create(0,1,eye_csr,MPI_COMM_SELF);
  las::Mat * col_mat = mat_fct->create(0,1,col_csr,MPI_COMM_SELF);
  assert(eye_mat && col_mat);
  auto * ops = las::getLASOps<las::sparse>();
  assert(ops);
  int rwcls[] = {0, 1, 2};
  ops->set(eye_mat,3,&rwcls[0],3,&rwcls[0],&eye[0]);
  ops->set(col_mat,3,&rwcls[0],3,&rwcls[0],&col[0]);
  las::MatMatMult * mlt = las::getSparseMatMatMult();
  las::Mat * rst = nullptr;
  mlt->exec(eye_mat,col_mat,&rst);
  double * rst_arr = nullptr;
  ops->get(rst,3,&rwcls[0],3,&rwcls[0],&rst_arr);
  for(int idx = 0; idx < 9; ++idx)
    assert(rst_arr[idx] == col[idx]);
  delete [] rst_arr;
  return 0;
}
