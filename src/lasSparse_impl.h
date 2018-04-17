#ifndef LAS_SPARSE_IMPL_H_
#define LAS_SPARSE_IMPL_H_
#include <cassert>
#include <cmath>
#include <cstring> //memset
#include "lasAlloc.h"
#include "lasCSR.h"
namespace las
{
  class csrMat
  {
    scalar * vls;
    CSR * csr;
  public:
    csrMat(CSR * c)
      : vls(nullptr)
      , csr(c)
    {
      alloc<Malloc>((void**)&vls,sizeof(scalar) * (c->getNumNonzero() + 1));
      memset(&vls[0],0,sizeof(scalar)*(csr->getNumNonzero()+1));
    }
    ~csrMat()
    {
      dealloc<Malloc>((void**)&vls);
    }
    scalar & operator()(int rr, int cc)
    {
      int idx[] = { (*csr)(rr,cc), csr->getNumNonzero() };
      bool dmy = rr < 0 || cc < 0;
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
      memset(&vls[0],0,sizeof(scalar)*(csr->getNumNonzero()+1));
    }
  };
  class simpleVec
  {
  private:
    scalar * vls;
    int cnt;
  public:
    simpleVec(int n)
      : vls(nullptr)
      , cnt(n)
    {
      alloc<Malloc>((void**)&vls,sizeof(scalar)*(n+1));
    }
    ~simpleVec()
    {
      dealloc<Malloc>((void**)vls);
    }
    scalar & operator[](int idx)
    {
      assert(idx < cnt);
      if(idx < 0)
        idx = cnt;
      return vls[idx];
    }
    int size()
    {
      return cnt;
    }
    // too coupled to the implementation to leave external
    void zero()
    {
      memset(&vls[0],0,sizeof(scalar)*(cnt+1));
    }
  };
  inline csrMat * getCSRMat(Mat * m)
  {
    return reinterpret_cast<csrMat*>(m);
  }
  inline simpleVec * getSimpleVec(Vec * v)
  {
    return reinterpret_cast<simpleVec*>(v);
  }
  inline Mat * createCSRMatrix(Sparsity * csr)
  {
    return reinterpret_cast<Mat*>(new csrMat(reinterpret_cast<CSR*>(csr)));
  }
  inline void destroyCSRMatrix(Mat * m)
  {
    delete getCSRMat(m);
  }
  inline Vec * createVector(unsigned n)
  {
    return reinterpret_cast<Vec*>(new simpleVec(n));
  }
  inline void destroyVector(Vec * v)
  {
    delete getSimpleVec(v);
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
  LasCreateMat * getMatBuilder<sparse>(int)
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
  LasCreateVec * getVecBuilder<sparse>(int)
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
      getSimpleVec(v)->zero();
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
      simpleVec * vec = getSimpleVec(v);
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
          if(vl != 0.0) // don't want to attempt to access zero locations in a sparse matrix
            (*mat)(rws[ii],cls[jj]) += vls[ii * cntc + jj];
        }
    }
    void _set(Vec * v, int cnt, int * rws, scalar * vls)
    {
      simpleVec * vec = getSimpleVec(v);
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
      simpleVec * vec = getSimpleVec(v);
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
          (*vls)[ii*cntc+jj] = (*mat)(rws[ii],cls[ii]);
    }
    scalar _norm(Vec * v)
    {
      simpleVec * vec = getSimpleVec(v);
      scalar nrm = 0.0;
      for(int ii = 0; ii < vec->size(); ++ii)
        nrm += (*vec)[ii] * (*vec)[ii];
      nrm = sqrt(nrm);
      return nrm;
    }
    scalar _dot(Vec * v0, Vec * v1)
    {
      simpleVec * vec0 = getSimpleVec(v0);
      simpleVec * vec1 = getSimpleVec(v1);
      int sz0 = vec0->size();
      assert(sz0 == vec1->size());
      scalar dt = 0.0;
      for(int ii = 0; ii < sz0; ++ii)
        dt += (*vec0)[ii] * (*vec1)[ii];
      return dt;
    }
    void _axpy(scalar a, Vec * x, Vec * y)
    {
      simpleVec * vx = getSimpleVec(x);
      simpleVec * vy = getSimpleVec(y);
      int szx = vx->size();
      assert(szx == vy->size());
      for(int ii = 0; ii < szx; ++ii)
        (*vy)[ii] = a * (*vx)[ii] + (*vy)[ii];

    }
    void _get(Vec * v, scalar *& vls)
    {
      simpleVec * vec = getSimpleVec(v);
      vls = &(*vec)[0];
    }
    void _restore(Vec*, scalar *&)
    { }
  };
  template <>
  inline LasOps<sparse> * getLASOps()
  {
    static sparse * ops = nullptr;
    if(ops == nullptr)
      ops = new sparse;
    return ops;
  }
}
#endif
