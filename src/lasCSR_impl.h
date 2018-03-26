#ifndef LAS_CSR_IMPL_H_
#define LAS_CSR_IMPL_H_
namespace las
{
  inline CSR::CSR(int ne, int nz, int * rs, int * cs)
    : neq(ne)
    , nnz(nz)
    , rws(rs,rs+ne+1)
    , cls(cs,cs+nz)
  { }
  /*
  inline int CSR::operator()(int rw, int cl) const
  {
  }
  */
}
#endif
