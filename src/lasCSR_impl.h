#ifndef LAS_CSR_IMPL_H_
#define LAS_CSR_IMPL_H_
namespace las
{
  inline CSR::CSR(int r, int c,  int nz, int * rs, int * cs)
    : nr(r)
    , nc(c)
    , nnz(nz)
    , rws(rs,rs+nr+1)
    , cls(cs,cs+nz)
  { }
}
#endif
