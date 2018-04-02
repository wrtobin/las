#include "lasCSR.h"
#include <cstring>
#include <iostream>
namespace las
{
  void constructFullMatrix(CSR * csr, scalar * sprs_mat, scalar * fll_mat)
  {
    int neq = csr->getNumEqs();
    memset(&fll_mat[0],0,neq*neq*sizeof(scalar));
    int lc = -1;
    for(int ii = 0; ii < neq; ii++)
      for(int jj = 0; jj < neq; jj++)
        if((lc = (*csr)(ii,jj)) != -1)
          fll_mat[ii*neq + jj] = sprs_mat[lc];
  }
}
