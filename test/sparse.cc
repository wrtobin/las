#include "lasSparse.h"
#include <cassert>
#include <sstream>
#include <mpi.h>
#include <las.h>
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
  las::MatMatMult * mlt = las::getMatMatMult<las::sparse>();
  las::Mat * rst = nullptr;
  mlt->exec(eye_mat,col_mat,&rst);
  double * rst_arr = nullptr;
  ops->get(rst,3,&rwcls[0],3,&rwcls[0],&rst_arr);
  for(int idx = 0; idx < 9; ++idx)
    assert(rst_arr[idx] == col[idx]);
  delete [] rst_arr;
  mat_fct->destroy(rst);
  double empty_row[] = {1.0, 0.0, 1.0,
                        0.0, 0.0, 0.0,
                        1.0, 0.0, 1.0};
  double row_col[] = {1.0, 1.0, 1.0,
                      1.0, 0.0, 0.0,
                      1.0, 0.0, 0.0};
  las::Sparsity * empty_row_csr = las::csrFromFull(&empty_row[0],3,3);
  las::Sparsity * row_col_csr = las::csrFromFull(&row_col[0],3,3);
  las::Mat * empty_row_mat = mat_fct->create(LAS_IGNORE,LAS_IGNORE,empty_row_csr,MPI_COMM_SELF);
  las::Mat * row_col_mat = mat_fct->create(LAS_IGNORE,LAS_IGNORE,row_col_csr,MPI_COMM_SELF);
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
  las::Sparsity * scalar_mat_mult_csr = las::csrFromFull(&test_scalarMatMult_vls[0], 3, 3);
  las::Mat * test_scalarMatMult = mat_fct->create(
      LAS_IGNORE, LAS_IGNORE,
      scalar_mat_mult_csr, MPI_COMM_SELF);
  ops->set(test_scalarMatMult, 3, &rwcls[0], 3, &rwcls[0],
           &test_scalarMatMult_vls[0]);
  las::Mat * scalarMatMult_rslt = nullptr;
  las::ScalarMatMult * smm = las::getScalarMatMult<las::sparse>();
  smm->exec(2.0, test_scalarMatMult, &scalarMatMult_rslt);
  double * scalarMatMult_rslt_arr = nullptr;
  ops->get(scalarMatMult_rslt, 3, &rwcls[0], 3, &rwcls[0],
           &scalarMatMult_rslt_arr);
  for (int i = 0; i < 9; ++i)
  {
    assert(fabs(test_scalarMatMult_vls[i] * 2 - scalarMatMult_rslt_arr[i]) <
           1E-15);
  }
  delete [] scalarMatMult_rslt_arr;
  scalarMatMult_rslt_arr = nullptr;
  mat_fct->destroy(scalarMatMult_rslt);
  // test inplace scalar mat mult
  smm->exec(2.0, test_scalarMatMult, nullptr);
  ops->get(test_scalarMatMult, 3, &rwcls[0], 3, &rwcls[0],
           &scalarMatMult_rslt_arr);
  for (int i = 0; i < 9; ++i)
  {
    assert(fabs(test_scalarMatMult_vls[i] * 2 - scalarMatMult_rslt_arr[i]) <
           1E-15);
  }
  mat_fct->destroy(test_scalarMatMult);
  delete [] scalarMatMult_rslt_arr;
  scalarMatMult_rslt_arr = nullptr;
  las::MatMatAdd * smsma = las::getMatMatAdd<las::sparse>();
  double smsma_expected_1[9] = {1.0 ,3.14, 0.0, 0.0, 4.14, 0.0, 0.0, 3.14,1.0};
  double smsma_expected_2[9] = {2.0 ,9.42, 0.0, 0.0, 11.42, 0.0, 0.0, 9.42,2.0};
  las::Mat * c=0;
  smsma->exec(1, eye_mat, 1, col_mat, &c);
  for(int i=0; i<3; ++i)
  {
    for(int j=0; j<3; ++j)
    {
    assert(fabs(smsma_expected_1[i*3+j] - (*las::getCSRMat(c))(i,j)) < 1E-15);
    }
  }
  smsma->exec(2.0, eye_mat, 3.0, col_mat, &c);
  for(int i=0; i<3; ++i)
  {
    for(int j=0; j<3; ++j)
    {
    assert(fabs(smsma_expected_2[i*3+j] - (*las::getCSRMat(c))(i,j)) < 1E-15);
    }
  }
  las::LasCreateVec * vb = las::getVecBuilder<las::sparse>(0);
  las::Vec* v1 =  vb->createRHS(eye_mat);
  las::Vec* v2 =  vb->createRHS(eye_mat);
  double v1_vals[3] = {1.0, 2.1, 3.0};
  double v2_vals[3] = {4.2, 5.0, 6.0};
  ops->set(v1, 3, &rwcls[0], &v1_vals[0]);
  ops->set(v2, 3, &rwcls[0], &v2_vals[0]);
  las::VecVecAdd * vva = las::getVecVecAdd<las::sparse>();
  las::Vec * v3 = nullptr;
  double * vva_rslts = nullptr;
  vva->exec(1, v1, 1, v2, v3);
  ops->get(v3, vva_rslts);
  for(int i=0; i<3; ++i)
  {
    assert(fabs(vva_rslts[i]-(1.0*v1_vals[i]+1.0*v2_vals[i])) < 1E-15);
  }
  vva->exec(0.5, v1, 2.0, v2, v3);
  for(int i=0; i<3; ++i)
  {
    assert(fabs(vva_rslts[i]-(0.5*v1_vals[i]+2.0*v2_vals[i])) < 1E-15);
  }
  vb->destroy(v1);
  vb->destroy(v2);
  vb->destroy(v3);
  las::destroySparsity<las::sparse>(eye_csr); 
  las::destroySparsity<las::sparse>(col_csr); 
  las::destroySparsity<las::sparse>(empty_row_csr); 
  las::destroySparsity<las::sparse>(row_col_csr); 
  las::destroySparsity<las::sparse>(scalar_mat_mult_csr); 
  mat_fct->destroy(c);
  mat_fct->destroy(eye_mat);
  mat_fct->destroy(col_mat);
  MPI_Finalize();
  return 0;
}
