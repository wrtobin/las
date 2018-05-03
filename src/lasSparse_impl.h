#ifndef LAS_SPARSE_IMPL_H_
#define LAS_SPARSE_IMPL_H_
#include "lasAlloc.h"
#include "lasCSR.h"
#include "lasInline.h"
#include "lasVec_impl.h"
#include <cassert>
#include <cmath>
#include <cstring> //memset
namespace las
{
  class csrMat
  {
    scalar * vls; // nnz +2 (dumy write, 0.0 read)
    CSR * csr;
  public:
    csrMat(CSR * c)
      : vls(nullptr)
      , csr(c)
    {
      alloc<Malloc>((void**)&vls,sizeof(scalar) * (c->getNumNonzero() + 2));
      memset(&vls[0],0,sizeof(scalar)*(csr->getNumNonzero() + 2));
    }
    ~csrMat()
    {
      dealloc<Malloc>((void**)&vls);
    }
    scalar & operator()(int rr, int cc)
    {
      int idx[] = { (*csr)(rr,cc), csr->getNumNonzero(), csr->getNumNonzero()+1 };
      int dmy = (rr < 0 || cc < 0) ? 1 : idx[0] < 0 ? 2 : 0;
      return vls[idx[dmy]];
    }
    CSR * getCSR()
    {
      return csr;
    }
    scalar * getVals()
    {
      return &vls[0];
    }
    // too coupled to the implementation to leave external
    void zero()
    {
      memset(&vls[0],0,sizeof(scalar)*(csr->getNumNonzero() + 2));
    }
  };
  LAS_INLINE csrMat * getCSRMat(Mat * m)
  {
    return reinterpret_cast<csrMat*>(m);
  }
  LAS_INLINE Mat * createCSRMatrix(Sparsity * csr)
  {
    return reinterpret_cast<Mat*>(new csrMat(reinterpret_cast<CSR*>(csr)));
  }
  LAS_INLINE void destroyCSRMatrix(Mat * m)
  {
    delete getCSRMat(m);
  }
  class csrMatBuilder : public LasCreateMat
  {
  public:
    virtual ~csrMatBuilder() {};
    Mat * create(unsigned,unsigned,Sparsity * s, MPI_Comm)
    {
      return createCSRMatrix(s);
    }
    void destroy(Mat * m)
    {
      destroyCSRMatrix(m);
    }
  };
  template <>
  LAS_INLINE LasCreateMat * getMatBuilder<sparse>(int)
  {
    static csrMatBuilder * mb = nullptr;
    if(mb == nullptr)
      mb = new csrMatBuilder;
    return mb;
  }
  class csrVecBuilder : public LasCreateVec
  {
  public:
    virtual ~csrVecBuilder() {};
    virtual Vec * create(unsigned lcl,unsigned,MPI_Comm)
    {
      return createVector(lcl);
    }
    virtual void destroy(Vec * v)
    {
      destroyVector(v);
    }
    virtual Vec * createLHS(Mat * m)
    {
      csrMat * cm = getCSRMat(m);
      CSR * csr = cm->getCSR();
      int cols = csr->getNumCols();
      return createVector(cols);
    }
    virtual Vec * createRHS(Mat * m)
    {
      csrMat * cm = getCSRMat(m);
      CSR * csr = cm->getCSR();
      int rows = csr->getNumRows();
      return createVector(rows);
    }
  };
  template <>
  LAS_INLINE LasCreateVec * getVecBuilder<sparse>(int)
  {
    static csrVecBuilder * vb = nullptr;
    if(vb == nullptr)
      vb = new csrVecBuilder;
    return vb;
  }
  class sparse : public LasOps<sparse>
  {
  public:
    void _zero(Mat * m)
    {
      getCSRMat(m)->zero();
    }
    void _zero(Vec * v)
    {
      getLASVec(v)->zero();
    }
    void _zero(Mat * m, int rw)
    {
      csrMat * mat = getCSRMat(m);
      int cols = mat->getCSR()->getNumCols();
      for(int ii = 0; ii < cols; ++ii)
        (*mat)(rw,ii) = 0.0;
    }
    void _assemble(Vec * v, int cnt, int * rws, scalar * vls)
    {
      lasVec * vec = getLASVec(v);
      for(int ii = 0; ii < cnt; ++ii)
        (*vec)[rws[ii]] += vls[ii];
    }
    void _assemble(Mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls)
    {
      csrMat * mat = getCSRMat(m);
      for(int ii = 0; ii < cntr; ++ii)
        for(int jj = 0; jj < cntc; ++jj)
        {
          scalar vl = vls[ii * cntc + jj];
          //if(vl != 0.0) // don't want to attempt to access zero locations in a sparse matrix
              (*mat)(rws[ii],cls[jj]) += vl;
        }
    }
    void _set(Vec * v, int cnt, int * rws, scalar * vls)
    {
      lasVec * vec = getLASVec(v);
      for(int ii = 0; ii < cnt; ++ii)
        (*vec)[rws[ii]] = vls[ii];
    }
    void _set(Mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls)
    {
      csrMat * mat = getCSRMat(m);
      for(int ii = 0; ii < cntr; ++ii)
        for(int jj = 0; jj < cntc; ++jj)
          (*mat)(rws[ii],cls[jj]) = vls[ii * cntc + jj];
    }
    void _get(Vec * v, int cntr, int * rws, scalar ** vls)
    {
      lasVec * vec = getLASVec(v);
      *vls = new scalar[cntr]();
      for(int ii = 0; ii < cntr; ++ii)
        (*vls)[ii] = (*vec)[rws[ii]];
    }
    void _get(Mat * m, int cntr, int * rws, int cntc, int * cls, scalar ** vls)
    {
      csrMat * mat = getCSRMat(m);
      *vls = new scalar[cntr*cntc]();
      for(int ii = 0; ii < cntr; ++ii)
        for(int jj = 0; jj < cntc; ++jj)
          (*vls)[ii*cntc+jj] = (*mat)(rws[ii],cls[jj]);
    }
    scalar _norm(Vec * v)
    {
      lasVec * vec = getLASVec(v);
      scalar nrm = 0.0;
      for(int ii = 0; ii < vec->size(); ++ii)
        nrm += (*vec)[ii] * (*vec)[ii];
      nrm = sqrt(nrm);
      return nrm;
    }
    scalar _dot(Vec * v0, Vec * v1)
    {
      lasVec * vec0 = getLASVec(v0);
      lasVec * vec1 = getLASVec(v1);
      int sz0 = vec0->size();
      assert(sz0 == vec1->size());
      scalar dt = 0.0;
      for(int ii = 0; ii < sz0; ++ii)
        dt += (*vec0)[ii] * (*vec1)[ii];
      return dt;
    }
    void _axpy(scalar a, Vec * x, Vec * y)
    {
      lasVec * vx = getLASVec(x);
      lasVec * vy = getLASVec(y);
      int szx = vx->size();
      assert(szx == vy->size());
      for(int ii = 0; ii < szx; ++ii)
        (*vy)[ii] = a * (*vx)[ii] + (*vy)[ii];

    }
    void _get(Vec * v, scalar *& vls)
    {
      lasVec * vec = getLASVec(v);
      vls = &(*vec)[0];
    }
    void _restore(Vec*, scalar *&)
    { }
  };
  template <>
  LAS_INLINE LasOps<sparse> * getLASOps()
  {
    static sparse * ops = nullptr;
    if(ops == nullptr)
      ops = new sparse;
    return ops;
  }
}
#endif
