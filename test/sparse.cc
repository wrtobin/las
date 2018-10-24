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
  // test scalarMat multiply
  double test_scalarMatMult_vls[] = {1.0, 2.0, 3.0, 4.0, 5.0,
                                     6.0, 7.0, 8.9, 9.0};
  las::Mat * test_scalarMatMult = mat_fct->create(
      LAS_IGNORE, LAS_IGNORE,
      las::csrFromFull(&test_scalarMatMult_vls[0], 3, 3), MPI_COMM_SELF);
  ops->set(test_scalarMatMult, 3, &rwcls[0], 3, &rwcls[0],
           &test_scalarMatMult_vls[0]);
  las::Mat * scalarMatMult_rslt = nullptr;
  las::ScalarMatMult * smm = las::getSparseScalarMatMult();
  smm->exec(2.0, test_scalarMatMult, &scalarMatMult_rslt);
  double * scalarMatMult_rslt_arr = nullptr;
  ops->get(scalarMatMult_rslt, 3, &rwcls[0], 3, &rwcls[0],
           &scalarMatMult_rslt_arr);
  for (int i = 0; i < 9; ++i)
  {
    assert(fabs(test_scalarMatMult_vls[i] * 2 - scalarMatMult_rslt_arr[i]) <
           1E-15);
  }
  // test inplace scalar mat mult
  smm->exec(2.0, test_scalarMatMult, nullptr);
  ops->get(test_scalarMatMult, 3, &rwcls[0], 3, &rwcls[0],
           &scalarMatMult_rslt_arr);
  for (int i = 0; i < 9; ++i)
  {
    assert(fabs(test_scalarMatMult_vls[i] * 2 - scalarMatMult_rslt_arr[i]) <
           1E-15);
  }
  MPI_Finalize();
  return 0;
}
