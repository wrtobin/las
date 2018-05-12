#include "lasSparse.h"
#include <cassert>
#include <sstream>
#include <mpi.h>
int main(int ac,char * av[])
{
  MPI_Init(&ac,&av);
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
  mat_fct->destroy(rst);
  mat_fct->destroy(eye_mat);
  mat_fct->destroy(col_mat);
  double empty_row[] = {1.0, 0.0, 1.0,
                        0.0, 0.0, 0.0,
                        1.0, 0.0, 1.0};
  double row_col[] = {1.0, 1.0, 1.0,
                      1.0, 0.0, 0.0,
                      1.0, 0.0, 0.0};
  las::Mat * empty_row_mat = mat_fct->create(LAS_IGNORE,LAS_IGNORE,las::csrFromFull(&empty_row[0],3,3),MPI_COMM_SELF);
  las::Mat * row_col_mat = mat_fct->create(LAS_IGNORE,LAS_IGNORE,las::csrFromFull(&row_col[0],3,3),MPI_COMM_SELF);
  ops->set(empty_row_mat,3,&rwcls[0],3,&rwcls[0],&empty_row[0]);
  ops->set(row_col_mat,3,&rwcls[0],3,&rwcls[0],&row_col[0]);
  las::Mat * rst2 = nullptr;
  mlt->exec(empty_row_mat,row_col_mat,&rst2);
  double * rst2_arr = nullptr;
  ops->get(rst2,3,&rwcls[0],3,&rwcls[0],&rst2_arr);
  double expect[] = {2.0, 1.0, 1.0,
                     0.0, 0.0, 0.0,
                     2.0, 1.0, 1.0};
  for(int idx = 0; idx < 9; idx++)
    assert(rst2_arr[idx] == expect[idx]);
  delete [] rst2_arr;
  mat_fct->destroy(rst2);
  mat_fct->destroy(row_col_mat);
  mat_fct->destroy(empty_row_mat);
  MPI_Finalize();
  return 0;
}
