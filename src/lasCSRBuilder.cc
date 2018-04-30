#include "lasCSRBuilder.h"
namespace las
{
  Sparsity * csrFromFull(scalar * mat, int rws, int cls)
  {
    CSRBuilder bld(rws,cls);
    for(int ii = 0; ii < rws; ii++)
      for(int jj = 0; jj < cls; jj++)
        if(mat[ii*cls + jj] != 0.0)
          bld.add(ii,jj);
    return bld.finalize();
  }
}
