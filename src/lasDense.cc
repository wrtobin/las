#include "lasDense.h"
namespace las
{
  class denseMatVec : public MatVecMult
  {
  public:
    void exec(Mat * x, Vec * a, Vec * b)
    {
      dnsMat * dns = getDnsMat(x);
      lasVec * va = getLASVec(a);
      lasVec * vb = getLASVec(b);
      int nr = dns->getNumRows();
      int nc = dns->getNumCols();
      int la = va->size();
      int lb = vb->size();
      assert(nc == la && "Matrix columns and lhs vector length must match");
      assert(nr == lb && "Matrix rows and rhs vector lenght must match");
      for(int rw = 0; rw < nr; ++rw)
        for(int cl = 0; cl < nc; ++cl)
          (*vb)[rw] += (*dns)(rw,cl) * (*va)[cl];
    }
  };
  class denseMatMat : public MatMatMult
  {
  public:
    void exec(Mat * xm, Mat * ym, Mat ** zm)
    {
      dnsMat * x = getDnsMat(xm);
      dnsMat * y = getDnsMat(ym);
      int x_nr = x->getNumRows();
      int x_nc = x->getNumCols();
      int y_nr = y->getNumRows();
      int y_nc = y->getNumCols();
      assert(x_nc == y_nr && "Matrix X cols must equal matrix Y rows");
      Sparsity * dnsty = createDensity(x_nr,y_nc);
      (*zm) = createDnsMat(dnsty);
      delete reinterpret_cast<Density*>(dnsty);
      dnsMat * z = getDnsMat(*zm);
      for(int xrw = 0; xrw < x_nr; ++xrw)
        for(int ycl = 0; ycl < y_nc; ++ycl)
          for(int inr = 0; inr < x_nc; ++inr)
            (*z)(xrw,ycl) += (*x)(xrw,inr) * (*y)(inr,ycl);
    }
  };
  MatVecMult * getDenseMatVecMult()
  {
    static denseMatVec * mv = nullptr;
    if(mv == nullptr)
      mv = new denseMatVec;
    return mv;
  }
  MatMatMult * getDenseMatMatMult()
  {
    static denseMatMat * mm = nullptr;
    if(mm == nullptr)
      mm = new denseMatMat;
    return mm;
  }
}
