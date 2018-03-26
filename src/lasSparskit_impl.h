#ifndef LAS_SPARSKIT_IMPL_H_
#define LAS_SPARSKIT_IMPL_H_
#include "lasSparskitExterns.h"
#include "lasCSR.h"
#include "lasDebug.h"
#include <cassert>
namespace las
{
  class skMat
  {
  private:
    double * vls;
    CSR * csr;
  public:
    skMat(CSR * c)
      : vls(new double [c->getNumNonzero()+1])
      , csr(c)
    {
      memset(&vls[0],0.0,sizeof(double)*(csr->getNumNonzero()+1));
    }
    ~skMat()
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
    // too coupled to the implementation to leave external
    void zero()
    {
      memset(&vls[0],0.0,sizeof(double)*(csr->getNumNonzero()+1));
    }
  };
  class skVec
  {
  private:
    double * vls;
    int cnt;
  public:
    skVec(int n)
      : vls(new double [n+1])
      , cnt(n)
    { }
    ~skVec()
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
  inline skMat * getSparskitMatrix(Mat * m)
  {
    return reinterpret_cast<skMat*>(m);
  }
  inline skVec * getSparskitVector(Vec * v)
  {
    return reinterpret_cast<skVec*>(v);
  }
  class SparskitLU : public LasSolve
  {
  protected:
    SparskitBuffers * bfrs;
    friend class SparskitQuickLU;
  public:
    SparskitLU(SparskitBuffers * b) : bfrs(b) {}
    virtual void solve(Mat * k, Vec * u, Vec * f);
  };
  // only perform the solve, do not decompose the matrix
  class SparskitQuickLU : public SparskitLU
  {
  public:
    SparskitQuickLU(SparskitBuffers * b) : SparskitLU(b) {}
    SparskitQuickLU(SparskitLU * lu) : SparskitLU(lu->bfrs) {}
    // the matrix k must have a csr format identical to that used previously in a normal SparskitLU solve
    virtual void solve(Mat * k, Vec * u, Vec * f);
  };
  inline void skOps::_zero(Mat * m)
  {
    getSparskitMatrix(m)->zero();
  }
  inline void skOps::_zero(Vec * v)
  {
    getSparskitVector(v)->zero();
  }
  inline void skOps::_assemble(Vec * v, int cnt, int * rws, double * vls)
  {
    skVec * vec = getSparskitVector(v);
    for(int ii = 0; ii < cnt; ++ii)
      (*vec)[rws[ii]] += vls[ii];
  }
  inline void skOps::_assemble(Mat * m, int cntr, int * rws, int cntc, int * cols, double * vls)
  {
    skMat * mat = getSparskitMatrix(m);
    for(int ii = 0; ii < cntr; ++ii)
      for(int jj = 0; jj < cntc; ++jj)
      {
        double vl = vls[ii * cntc + jj];
        if(vl != 0.0) // don't want to attempt to access zero locations in a sparse matrix
          (*mat)(rws[ii],cols[jj]) += vls[ii * cntc + jj];
      }
  }
  inline void skOps::_set(Vec * v, int cnt, int * rws, double * vls)
  {
    skVec * vec = getSparskitVector(v);
    for(int ii = 0; ii < cnt; ++ii)
      (*vec)[rws[ii]] = vls[ii];
  }
  inline void skOps::_set(Mat * m, int cntr, int * rws, int cntc, int * cols, double * vls)
  {
    skMat * mat = getSparskitMatrix(m);
    for(int ii = 0; ii < cntr; ++ii)
      for(int jj = 0; jj < cntc; ++jj)
        (*mat)(rws[ii],cols[jj]) = vls[ii * cntc + jj];
  }
  inline double skOps::_norm(Vec * v)
  {
    skVec * vec = getSparskitVector(v);
    double nrm = 0.0;
    for(int ii = 0; ii < vec->size(); ++ii)
      nrm += (*vec)[ii] * (*vec)[ii];
    nrm = sqrt(nrm);
    return nrm;
  }
  inline double skOps::_dot(Vec * v0, Vec * v1)
  {
    skVec * vec0 = getSparskitVector(v0);
    skVec * vec1 = getSparskitVector(v1);
    int sz0 = vec0->size();
    assert(sz0 == vec1->size());
    double dt = 0.0;
    for(int ii = 0; ii < sz0; ++ii)
      dt += (*vec0)[ii] * (*vec1)[ii];
    return dt;
  }
  // y = ax + y
  inline void skOps::_axpy(double a, Vec * x, Vec * y)
  {
    skVec * vx = getSparskitVector(x);
    skVec * vy = getSparskitVector(y);
    int szx = vx->size();
    assert(szx == vy->size());
    for(int ii = 0; ii < szx; ++ii)
      (*vy)[ii] = a * (*vx)[ii] + (*vy)[ii];
  }
  inline void skOps::_get(Vec * v, double *& vls)
  {
    skVec * vec = getSparskitVector(v);
    vls = &(*vec)[0];
  }
  inline Mat * createSparskitMatrix(CSR * csr)
  {
    return reinterpret_cast<Mat*>(new skMat(csr));
  }
  inline void deleteSparskitMatrix(Mat * m)
  {
    skMat * mat = getSparskitMatrix(m);
    delete mat;
  }
  inline Vec * createSparskitVector(int n)
  {
    return reinterpret_cast<Vec*>(new skVec(n));
  }
  inline void deleteSparskitVector(Vec * v)
  {
    skVec * vec = getSparskitVector(v);
    delete vec;
  }
  inline LasOps<skOps> * initSparskitOps()
  {
    static skOps * ops = NULL;
    if(ops == NULL)
      ops = new skOps;
    return ops;
  }
  inline LasSolve * createSparskitLUSolve(SparskitBuffers * b)
  {
    return new SparskitLU(b);
  }
  inline LasSolve * createSparskitQuickLUSolve(SparskitBuffers * b)
  {
    return new SparskitQuickLU(b);
  }
  inline LasSolve * createSparskitQuickLUSolve(LasSolve * slv)
  {
    SparskitLU * skt_slv = reinterpret_cast<SparskitLU*>(slv);
    return new SparskitQuickLU(skt_slv);
  }
  inline void printSparskitMat(std::ostream & o, Mat * mi)
  {
    skMat * m = getSparskitMatrix(mi);
    int ndofs = m->getCSR()->getNumEqs();
    for(int rr = 0; rr < ndofs; ++rr)
    {
      for(int cc = 0; cc < ndofs; ++cc)
      {
        o << (*m)(rr,cc) << ' ';
      }
      o << '\b' << std::endl;
    }
  }
  inline double getSparskitMatValue(Mat * k, int rr, int cc)
  {
    skMat * m = getSparskitMatrix(k);
    DBG(int ndofs = m->getCSR()->getNumEqs());
    assert(rr < ndofs && rr >= 0);
    assert(cc < ndofs && cc >= 0);
    return (*m)(rr,cc);
  }
  inline void setSparskitMatValue(Mat * k, int rr, int cc, double vl)
  {
    skMat * m = getSparskitMatrix(k);
    DBG(int ndofs = m->getCSR()->getNumEqs());
    assert(rr < ndofs && rr >= 0);
    assert(cc < ndofs && cc >= 0);
    (*m)(rr,cc) = vl;
  }
  inline void SparskitLU::solve(Mat * k, Vec * u, Vec * f)
  {
    bfrs->zero();
    double tol = 1e-6;
    skMat * mat = getSparskitMatrix(k);
    skVec * uv = getSparskitVector(u);
    skVec * fv = getSparskitVector(f);
    CSR * csr = mat->getCSR();
    int bfr_lng = bfrs->matrixLength();
    int ierr = 0;
    int ndofs = csr->getNumEqs();
    ilut_(&ndofs,
          &(*mat)(0,0),
          csr->getCols(),
          csr->getRows(),
          &ndofs,
          &tol,
          bfrs->matrixBuffer(),
          bfrs->colsBuffer(),
          bfrs->rowsBuffer(),
          &bfr_lng,
          bfrs->doubleWorkBuffer(),
          bfrs->intWorkBuffer(),
          &ierr);
    if(ierr != 0)
    {
      std::cerr << "ERROR: ilut_ returned error code " << ierr << std::endl;
      return;
    }
    lusol_(&ndofs,
           &(*fv)[0],
           &(*uv)[0],
           bfrs->matrixBuffer(),
           bfrs->colsBuffer(),
           bfrs->rowsBuffer());
  }
  inline void SparskitQuickLU::solve(Mat * k, Vec * u, Vec * f)
  {
    skMat * mat = getSparskitMatrix(k);
    skVec * uv = getSparskitVector(u);
    skVec * fv = getSparskitVector(f);
    int ndofs = mat->getCSR()->getNumEqs();
    lusol_(&ndofs,
           &(*fv)[0],
           &(*uv)[0],
           bfrs->matrixBuffer(),
           bfrs->colsBuffer(),
           bfrs->rowsBuffer());
  }
}
#endif
