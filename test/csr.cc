#include "lasSparse.h"
#include <cassert>
#include <sstream>
int main(int argc, char * argv[])
{
  double eye[] = { 1.0 , 0.0 , 0.0 ,
                   0.0 , 1.0 , 0.0 ,
                   0.0 , 0.0 , 1.0 };
  las::CSR * eye_csr = las::csrFromFull(&eye[0],3,3);
  assert(eye_csr);
  las::Mat * mat_csr = las::createCSRMatrix(eye_csr);
  assert(mat_csr);
  las::LasOps<las::csrOps> * ops = las::initCSROps();
  assert(ops);
  int rwcls[] = {0, 1, 2};
  double vl = 1.0;
  ops->set(mat_csr,1,&rwcls[0],1,&rwcls[0],&vl);
  ops->set(mat_csr,1,&rwcls[1],1,&rwcls[1],&vl);
  ops->set(mat_csr,1,&rwcls[2],1,&rwcls[2],&vl);
  return 0;
}
