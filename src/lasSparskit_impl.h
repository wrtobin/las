#ifndef LAS_SPARSKIT_IMPL_H_
#define LAS_SPARSKIT_IMPL_H_
#include "lasSparskitExterns.h"
#include "lasSparse_impl.h"
#include "lasDebug.h"
#include <cassert>
#include <iostream>
namespace las
{
  typedef csrMat skMat;
  typedef simpleVec skVec;
  inline skMat * getSparskitMatrix(Mat * m)
  {
    return reinterpret_cast<skMat*>(m);
  }
  inline skVec * getSparskitVector(Vec * v)
  {
    return reinterpret_cast<skVec*>(v);
  }
  class SparskitLU : public Solve
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
  inline Solve * createSparskitLUSolve(SparskitBuffers * b)
  {
    return new SparskitLU(b);
  }
  inline Solve * createSparskitQuickLUSolve(SparskitBuffers * b)
  {
    return new SparskitQuickLU(b);
  }
  inline Solve * createSparskitQuickLUSolve(Solve * slv)
  {
    SparskitLU * skt_slv = reinterpret_cast<SparskitLU*>(slv);
    return new SparskitQuickLU(skt_slv);
  }
  inline void printSparskitMat(std::ostream & o, Mat * mi)
  {
    skMat * m = getSparskitMatrix(mi);
    int ndofs = m->getCSR()->getNumRows();
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
    DBG(int ndofs = m->getCSR()->getNumRows());
    assert(rr < ndofs && rr >= 0);
    assert(cc < ndofs && cc >= 0);
    return (*m)(rr,cc);
  }
  inline void setSparskitMatValue(Mat * k, int rr, int cc, double vl)
  {
    skMat * m = getSparskitMatrix(k);
    DBG(int ndofs = m->getCSR()->getNumRows());
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
    int ndofs = csr->getNumRows();
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
    int ndofs = mat->getCSR()->getNumRows();
    lusol_(&ndofs,
           &(*fv)[0],
           &(*uv)[0],
           bfrs->matrixBuffer(),
           bfrs->colsBuffer(),
           bfrs->rowsBuffer());
  }
}
#endif
