#ifndef LAS_H_
#define LAS_H_
#include "lasComm.h"
#include "lasScalar.h"
namespace las
{
  class Mat;
  class Vec;
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
      static_cast<T*>(this)->_norm(v);
    }
    scalar dot(Vec * v0, Vec * v1)
    {
      static_cast<T*>(this)->_dot(v0,v1);
    }
    void axpy(scalar a, Vec * x, Vec * y)
    {
      static_cast<T*>(this)->_axpy(a,x,y);
    }
    void get(Vec * v, scalar *& vls)
    {
      static_cast<T*>(this)->_vet(v,vls);
    }
    void restore(Vec * v, scalar *& vls)
    {
      static_cast<T*>(this)->_restore(v,vls);
    }
  };
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
