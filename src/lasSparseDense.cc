#include "lasSparseDense.h"
namespace las
{
  class sparseMatDenseMat : public MatMatMult
  {
  public:
    void exec(Mat * xm, Mat * ym, Mat ** zm)
    {
      csrMat * x = getCSRMat(xm);
      dnsMat * y = getDnsMat(ym);
      CSR * x_csr = x->getCSR();
      int x_nr = x_csr->getNumRows();
      int x_nc = x_csr->getNumCols();
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
  MatMatMult * getSparseMatDenseMatMult()
  {
    static sparseMatDenseMat * mm = nullptr;
    if(mm == nullptr)
      mm = new sparseMatDenseMat;
    return mm;
  }
}
