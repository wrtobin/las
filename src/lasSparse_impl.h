#ifndef LAS_SPARSE_IMPL_H_
#define LAS_SPARSE_IMPL_H_
#include <cassert>
#include <cmath>
namespace las
{
  class csrMat
  {
    double * vls;
    CSR * csr;
  public:
    csrMat(CSR * c)
      : vls(new double [c->getNumNonzero()+1])
      , csr(c)
    {
      memset(&vls[0],0.0,sizeof(double)*(csr->getNumNonzero()+1));
    }
    ~csrMat()
    {
      delete [] vls;
    }
    double & operator()(int rr, int cc)
    {
      int idx = -1;
      if(rr < 0 || cc < 0)
        idx = csr->getNumNonzero();
      else
        idx = (*csr)(rr,cc);
      return vls[idx];
    }
    CSR * getCSR()
    {
      return csr;
    }
    double * getVals()
    {
      return &vls[0];
    }
    // too coupled to the implementation to leave external
    void zero()
    {
      memset(&vls[0],0.0,sizeof(double)*(csr->getNumNonzero()+1));
    }
  };
  class simpleVec
  {
  private:
    double * vls;
    int cnt;
  public:
    simpleVec(int n)
      : vls(new double [n+1])
      , cnt(n)
    { }
    ~simpleVec()
    {
      delete [] vls;
    }
    double & operator[](int idx)
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
      memset(&vls[0],0.0,sizeof(double)*(cnt+1));
    }
  };
  inline LasOps<csrOps> * initCSROps()
  {
    static csrOps * ops = nullptr;
    if(ops == nullptr)
      ops = new csrOps;
    return ops;
  }
  inline csrMat * getCSRMat(Mat * m)
  {
    return reinterpret_cast<csrMat*>(m);
  }
  inline simpleVec * getSimpleVec(Vec * v)
  {
    return reinterpret_cast<simpleVec*>(v);
  }
  inline Mat * createCSRMatrix(CSR * csr)
  {
    return reinterpret_cast<Mat*>(new csrMat(csr));
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
  inline void csrOps::_zero(Mat * m)
  {
    getCSRMat(m)->zero();
  }
  inline void csrOps::_zero(Vec * v)
  {
    getSimpleVec(v)->zero();
  }
  inline void csrOps::_assemble(Vec * v, int cnt, int * rws, double * vls)
  {
    simpleVec * vec = getSimpleVec(v);
    for(int ii = 0; ii < cnt; ++ii)
      (*vec)[rws[ii]] += vls[ii];
  }
  inline void csrOps::_assemble(Mat * m, int cntr, int * rws, int cntc, int * cols, double * vls)
  {
    csrMat * mat = getCSRMat(m);
    for(int ii = 0; ii < cntr; ++ii)
      for(int jj = 0; jj < cntc; ++jj)
      {
        double vl = vls[ii * cntc + jj];
        if(vl != 0.0) // don't want to attempt to access zero locations in a sparse matrix
          (*mat)(rws[ii],cols[jj]) += vls[ii * cntc + jj];
      }
  }
  inline void csrOps::_set(Vec * v, int cnt, int * rws, double * vls)
  {
    simpleVec * vec = getSimpleVec(v);
    for(int ii = 0; ii < cnt; ++ii)
      (*vec)[rws[ii]] = vls[ii];
  }
  inline void csrOps::_set(Mat * m, int cntr, int * rws, int cntc, int * cols, double * vls)
  {
    csrMat * mat = getCSRMat(m);
    for(int ii = 0; ii < cntr; ++ii)
      for(int jj = 0; jj < cntc; ++jj)
        (*mat)(rws[ii],cols[jj]) = vls[ii * cntc + jj];
  }
  inline double csrOps::_norm(Vec * v)
  {
    simpleVec * vec = getSimpleVec(v);
    double nrm = 0.0;
    for(int ii = 0; ii < vec->size(); ++ii)
      nrm += (*vec)[ii] * (*vec)[ii];
    nrm = sqrt(nrm);
    return nrm;
  }
  inline double csrOps::_dot(Vec * v0, Vec * v1)
  {
    simpleVec * vec0 = getSimpleVec(v0);
    simpleVec * vec1 = getSimpleVec(v1);
    int sz0 = vec0->size();
    assert(sz0 == vec1->size());
    double dt = 0.0;
    for(int ii = 0; ii < sz0; ++ii)
      dt += (*vec0)[ii] * (*vec1)[ii];
    return dt;
  }
  // y = ax + y
  inline void csrOps::_axpy(double a, Vec * x, Vec * y)
  {
    simpleVec * vx = getSimpleVec(x);
    simpleVec * vy = getSimpleVec(y);
    int szx = vx->size();
    assert(szx == vy->size());
    for(int ii = 0; ii < szx; ++ii)
      (*vy)[ii] = a * (*vx)[ii] + (*vy)[ii];
  }
  inline void csrOps::_get(Vec * v, double *& vls)
  {
    simpleVec * vec = getSimpleVec(v);
    vls = &(*vec)[0];
  }
}
#endif
