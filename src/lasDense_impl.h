#ifndef LAS_DENSE_IMPL_H_
#define LAS_DENSE_IMPL_H_
#include "lasVec_impl.h"
#include <cmath>
namespace las
{
  struct Density
  {
    int rws;
    int cls;
    scalar * vls;
    Density(int r, int c, scalar * v)
      : rws(r)
      , cls(c)
      , vls(v)
    { }
  };
  class dnsMat
  {
    int rws;
    int cls;
    scalar * vls;
    bool own;
  public:
    dnsMat(int rows, int cols, scalar * vals = nullptr)
      : rws(rows)
      , cls(cols)
      , vls(vals)
      , own(false)
    {
      if(vls == nullptr)
      {
        vls = new scalar[rws*cls];
        own = true;
      }
    }
    ~dnsMat()
    {
      if(own) delete [] vls;
    }
    scalar & operator()(int rr, int cc)
    {
      return vls[rr*cls + cc];
    }
    int getNumRows() { return rws; }
    int getNumCols() { return cls; }
    scalar * getVals()
    {
      return &vls[0];
    }
    void zero()
    {
      memset(&vls[0],0,sizeof(scalar)*rws*cls);
    }
  };
  LAS_INLINE Sparsity * createDensity(int rws, int cls, scalar * vls)
  {
    return reinterpret_cast<Sparsity*>(new Density(rws,cls,vls));
  }
  LAS_INLINE dnsMat * getDnsMat(Mat * m)
  {
    return reinterpret_cast<dnsMat*>(m);
  }
  LAS_INLINE Mat * createDnsMat(Sparsity * s)
  {
    Density * dns = reinterpret_cast<Density*>(s);
    return reinterpret_cast<Mat*>(new dnsMat(dns->rws,dns->cls,dns->vls));
  }
  LAS_INLINE void destroyDnsMat(Mat * m)
  {
    delete getDnsMat(m);
  }
  class dnsMatBuilder : public LasCreateMat
  {
  public:
    virtual ~dnsMatBuilder() {};
    Mat * create(unsigned,unsigned,Sparsity * s, MPI_Comm)
    {
      return createDnsMat(s);
    }
    void destroy(Mat * m)
    {
      destroyDnsMat(m);
    }
  };
  template <>
  LAS_INLINE LasCreateMat * getMatBuilder<dense>(int)
  {
    static dnsMatBuilder * mb = nullptr;
    if(mb == nullptr)
      mb = new dnsMatBuilder;
    return mb;
  }
  class dnsVecBuilder : public LasCreateVec
  {
  public:
    virtual ~dnsVecBuilder() {};
    virtual Vec * create(unsigned lcl,unsigned,MPI_Comm)
    {
      return createVector(lcl);
    }
    virtual Vec * create(scalar * data, unsigned lcl, unsigned, MPI_Comm)
    {
      return createVector(data, lcl);
    }
    virtual void destroy(Vec * v)
    {
      destroyVector(v);
    }
    virtual Vec * createLHS(Mat * m)
    {
      dnsMat * dns = getDnsMat(m);
      int cols = dns->getNumCols();
      return createVector(cols);
    }
    virtual Vec * createRHS(Mat * m)
    {
      dnsMat * dns = getDnsMat(m);
      int rows = dns->getNumRows();
      return createVector(rows);
    }
  };
  template <>
  LAS_INLINE LasCreateVec * getVecBuilder<dense>(int)
  {
    static dnsVecBuilder * vb = nullptr;
    if(vb == nullptr)
      vb = new dnsVecBuilder;
    return vb;
  }
  class dense : public LasOps<dense>
  {
  public:
    void _zero(Mat * m)
    {
      getDnsMat(m)->zero();
    }
    void _zero(Vec * v)
    {
      getLASVec(v)->zero();
    }
    void _zero(Mat * m, int rw)
    {
      dnsMat * dns = getDnsMat(m);
      int cols = dns->getNumCols();
      for(int ii = 0; ii < cols; ++ii)
        (*dns)(rw,ii) = 0.0;
    }
    void _assemble(Vec * v, int cnt, int * rws, scalar * vls)
    {
      lasVec * vec = getLASVec(v);
      for(int ii = 0; ii < cnt; ++ii)
        (*vec)[rws[ii]] += vls[ii];
    }
    void _assemble(Mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls)
    {
      dnsMat * dns = getDnsMat(m);
      for(int ii = 0; ii < cntr; ++ii)
        for(int jj = 0; jj < cntc; ++jj)
          (*dns)(rws[ii],cls[jj]) += vls[ii*cntc + jj];
    }
    void _set(Vec * v, int cnt, int * rws, scalar * vls)
    {
      lasVec * vec = getLASVec(v);
      for(int ii = 0; ii < cnt; ++ii)
        (*vec)[rws[ii]] = vls[ii];
    }
    void _set(Mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls)
    {
      dnsMat * dns = getDnsMat(m);
      for(int ii = 0; ii < cntr; ++ii)
        for(int jj = 0; jj < cntc; ++jj)
          (*dns)(rws[ii],cls[jj]) = vls[ii*cntc + jj];
    }
    void _get(Vec * v, int cntr, int * rws, scalar ** vls)
    {
      lasVec * vec = getLASVec(v);
      *vls = new scalar[cntr]();
      for(int ii = 0; ii < cntr; ++ii)
        (*vls)[ii] = (*vec)[rws[ii]];
    }
    void _get(Mat * m, int cntr, int * rws, int  cntc, int * cls, scalar ** vls)
    {
      dnsMat * dns = getDnsMat(m);
      *vls = new scalar[cntr*cntc]();
      for(int ii = 0; ii < cntr; ++ii)
        for(int jj = 0; jj < cntc; ++jj)
          (*vls)[ii*cntc+jj] = (*dns)(rws[ii],cls[jj]);
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
      lasVec * a = getLASVec(v0);
      lasVec * b = getLASVec(v1);
      int sza = a->size();
      assert(sza == b->size());
      scalar dt = 0.0;
      for(int ii = 0; ii < sza; ++ii)
        dt += (*a)[ii] * (*b)[ii];
      return dt;
    }
    void _axpy(scalar a, Vec * vx, Vec * vy)
    {
      lasVec * x = getLASVec(vx);
      lasVec * y = getLASVec(vy);
      int szx = x->size();
      assert(szx == y->size());
      for(int ii = 0; ii < szx; ++ii)
        (*y)[ii] += a * (*x)[ii];
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
  LAS_INLINE LasOps<dense> * getLASOps()
  {
    static dense * ops = nullptr;
    if(ops == nullptr)
      ops = new dense;
    return ops;
  }
  
}
#endif
