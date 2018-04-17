#include "lasCSR.h"
#include <cassert>
namespace las
{
  bool needNonzero(int * rws, int rw, int * cls, int cl)
  {
    bool result = true;
    // negative rw and cl correspond to fixed dofs and won't add a nonzero
    if(rw < 0 || cl < 0)
      result = false;
    else
    {
      for(int ii = rws[rw]; ii < rws[rw+1]; ii++)
      {
        if(cls[ii-1] == (cl+1))
        {
          result = false;
          break;
        }
      }
    }
    return result;
  }
  int findLocation(int * rws, int rw, int * cls, int cl, int neq)
  {
    (void)neq;
    int loc = 0;
    assert(rw < neq);
    for(loc = rws[rw]; (loc < rws[rw+1]) && (cls[loc-1] < cl+1); ++loc) {}
    return loc;
  }
  bool addNonzero(int * rws, int rw, int * cls, int cl, int neq)
  {
    bool result = false;
    if((result = needNonzero(rws,rw,cls,cl)))
    {
      int ii = findLocation(rws,rw,cls,cl,neq) - 1;
      for(int jj = rw+1; jj < neq+1; ++jj)
        rws[jj]++;
      for(int jj = rws[neq] - 2; jj > ii; --jj)
        cls[jj] = cls[jj-1];
      cls[ii] = cl+1;
    }
    return result;
  }
  Sparsity * csrFromFull(scalar * mat, int rws, int cls)
  {
    int nnz = 0;
    std::vector<int> rwb(rws+1,1);
    std::vector<int> clb(rws*rws,0);
    for(int ii = 0; ii < rws; ii++)
      for(int jj = 0; jj < cls; jj++)
        if(mat[ii*cls + jj] != 0.0 && addNonzero(&rwb[0],ii,&clb[0],jj,rws))
          nnz++;
    rwb[rws+1] = nnz;
    CSR * rslt = new CSR(rws,cls,nnz,&rwb[0],&clb[0]);
    return reinterpret_cast<Sparsity*>(rslt);
  }
}
