#include "lasDense.h"
#include <cassert>
#include <vector>
#include <mpi.h>
int main(int ac, char * av[])
{
  MPI_Init(&ac,&av);
  int nr = 4;
  int nc = 3;
  double darr[] = {1.0, 2.0, 3.0,
                   1.0, 1.0, 1.0,
                   2.0, 2.0, 2.0,
                   3.0, 3.0, 3.0,
                   4.0, 4.0, 4.0,};
  las::LasCreateMat * dns_fct = las::getMatBuilder<las::dense>(0);
  las::Mat * dns_mat = dns_fct->create(0, 0, las::createDensity(nr, nc, darr), MPI_COMM_SELF);
  las::LasOps<las::dense> * dns_ops = las::getLASOps<las::dense>();
  std::vector <int> rws(nc, 0);
  std::vector <int> cls(nr, 0);
  for(int rw = 0; rw < nr; ++rw)
    rws[rw] = rw;
  for(int cl = 0; cl < nc; ++cl)
    cls[cl] = cl;
  double * val = nullptr;
  dns_ops->get(dns_mat, nr, &rws[0], nc, &cls[0], &val);
  for(int rw = 0; rw < nr; ++rw)
    for(int cl = 0; cl < nc; ++cl)
      assert(val[rw*nc + cl] == darr[rw*nc + cl]);
  delete [] val;
  MPI_Finalize();
  return 0;
}
