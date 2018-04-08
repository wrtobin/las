#ifndef LAS_H_
#define LAS_H_
#include "lasComm.h"
#include "lasScalar.h"
namespace las
{
  class Mat;
  class Vec;
  class Sparsity;
  template <class T>
  class LasOps
  {
  public:
    void zero(Mat * m)
    {
      static_cast<T*>(this)->_zero(m);
    }
    void zero(Vec * v)
    {
      static_cast<T*>(this)->_zero(v);
    }
    void zero(Mat * v, int rw)
    {
      static_cast<T*>(this)->_zero(v,rw);
    }
    void assemble(Vec * v, int cnt, int * rws, scalar * vls)
    {
      static_cast<T*>(this)->_assemble(v,cnt,rws,vls);
    }
    void assemble(Mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls)
    {
      static_cast<T*>(this)->_assemble(m,cntr,rws,cntc,cls,vls);
    }
    void set(Vec * v, int cnt, int * rws, scalar * vls)
    {
      static_cast<T*>(this)->_set(v,cnt,rws,vls);
    }
    void set(Mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls)
    {
      static_cast<T*>(this)->_set(m,cntr,rws,cntc,cls,vls);
    }
    scalar norm(Vec * v)
    {
      return static_cast<T*>(this)->_norm(v);
    }
    scalar dot(Vec * v0, Vec * v1)
    {
      return static_cast<T*>(this)->_dot(v0,v1);
    }
    void axpy(scalar a, Vec * x, Vec * y)
    {
      static_cast<T*>(this)->_axpy(a,x,y);
    }
    void get(Vec * v, scalar *& vls)
    {
      static_cast<T*>(this)->_get(v,vls);
    }
    void restore(Vec * v, scalar *& vls)
    {
      static_cast<T*>(this)->_restore(v,vls);
    }
  };
  class LasCreateMat
  {
  public:
    virtual Mat * create(unsigned lcl, unsigned bs, Sparsity * s, MPI_Comm cm) = 0;
    virtual void destroy(Mat * m) = 0;
  };
  class LasCreateVec
  {
  public:
    virtual Vec * create(unsigned lcl, unsigned bs, MPI_Comm cm) = 0;
    virtual void destroy(Vec * v) = 0;
    virtual Vec * createRHS(Mat * m);
    virtual Vec * createLHS(Mat * m);
  };
  template <typename T>
  LasOps<T> * getLASOps();
  template <typename T>
  LasCreateMat * getMatBuilder(int id);
  template <typename T>
  LasCreateVec * getVecBuilder(int id);
  class LasSolve
  {
  public:
    virtual void solve(Mat * k, Vec * u, Vec * f) = 0;
  };
  class LasMultiply
  {
  public:
    virtual void exec(Mat * x, Vec * a, Vec * b) = 0;
  };
}
#endif
