#include <las.h>
@backend_header_includes@
#include <cassert>
#include <iostream>
#define BACKEND @BACKEND@
bool close(double a, double b, double eps = 1e-8)
{
  return fabs(a-b) < fabs(eps);
}
int main(int argc, char * argv[])
{
  auto ops = las::getLASOps<BACKEND>();
  // test vector ops
  auto vb = las::getVecBuilder<BACKEND>();
  las::Vec * v0 = vb->create(3,1,LAS_COMM_WORLD);
  las::Vec * v1 = vb->create(3,1,LAS_COMM_WORLD);
  ops->zero(v0);
  double vls = {1.0, 1.0, 1.0 };
  int all_rws[] = {0, 1, 2};
  int evn_rws[] = {0, 2};
  las::scalar * v0_dat = nullptr;
  ops->get(v0,3,&all_rws[0],v0_dat);
  for(int rw = 0; rw < 3; ++rw)
  { assert(v0_dat[rw] == 0.0); }
  delete [] v0_dat;
  ops->assemble(v0,&all_rws[0],&vls[0]); // (111)
  ops->assemble(v0,&evn_rws[0],&vls[0]); // (212)
  ops->assemble(v0,1,&all_rws[2],&vls[0]); // (213)
  ops->get(v0,3,&all_rws[0],v0_dat);
  double rst = { 2.0, 1.0, 3.0 };
  for(int rw = 0; rw < 3; ++rw)
  { assert(close(v0_dat[rw],rst[rw])); }
  delete [] v0_dat;
  ops->set(v0,&all_rws[0],&vls[0]);
  ops->get(v0,3,&all_rws[0],v0_dat);
  for(int rw = 0; rw < 3; ++rw)
  { assert(close(v0_vat[rw],vls[rw])); }
  delete [] v0_dat;
  double vl = ops->norm(v0);
  assert(close(vl,sqrt(3)));
  ops->set(v1,&all_rws[0],&vls[0]);
  ops->assemble(v1,&all_rws[0],&vls[0]);
  vl = ops->dot(v0,v1);
  assert(close(vl,6.0));
  ops->axpy(2.0,v1,v0);
  ops->get(v0,3,&all_rws[0],v0_dat);
  double rst2 = { 5.0, 5.0, 5.0 };
  for(int rw = 0; rw < 3; ++rw)
  { assert(close(v0_dat[rw],rst2[rw])); }
  delete [] v0_dat;
  auto mb = las::getMatBuilder<BACKEND>();
}
