#include "lasSparse.h"
#include "lasCSRBuilder.h"
namespace las
{
  class sparseMatVec : public MatVecMult
  {
  public:
    void exec(Mat * x, Vec * a, Vec * b)
    {
      csrMat * cm = getCSRMat(x);
      CSR * csr = cm->getCSR();
      simpleVec * sa = getSimpleVec(a);
      simpleVec * sb = getSimpleVec(b);
      int nr = csr->getNumRows();
      int nc = csr->getNumCols();
      int la = sa->size();
      int lb = sb->size();
      assert(nc == la && "Matrix columns and lhs vector length must match");
      assert(nr == lb && "Matrix rows and rhs vector length must match");
      for(int rw = 0; rw < nr; ++rw)
        for(int cl = 0; cl < nc; ++cl)
          (*sb)[rw] += (*cm)(rw,cl) * (*sa)[cl];
    }
  };
  class sparseMatMat : public MatMatMult
  {
  public:
    void exec(Mat * x, Mat * y, Mat ** z)
    {
      csrMat * cx = getCSRMat(x);
      csrMat * cy = getCSRMat(y);
      CSR * x_csr = cx->getCSR();
      CSR * y_csr = cy->getCSR();
      int x_nr = x_csr->getNumRows();
      int x_nc = x_csr->getNumCols();
      int y_nr = y_csr->getNumRows();
      int y_nc = y_csr->getNumCols();
      assert(x_nc == y_nr && "Matrix X cols must equal Matrix Y rows");
      std::vector<int> rows(x_nr + 1,1);
      std::vector<int> cols(x_nr * y_nc,0);
      int nnz = 0;
      for(int xrw = 0; xrw < x_nr; ++xrw)
        for(int ycl = 0; ycl < y_nc; ++ycl)
          for(int inr = 0; inr < x_nc; ++inr)
            if(((*x_csr)(xrw,inr) >= 0 && (*y_csr)(inr,ycl) >= 0))
              if(addNonzero(&rows[0],xrw,&cols[0],ycl,x_nr))
                nnz++;
      CSR * z_csr = new CSR(x_nr,y_nc,nnz,&rows[0],&cols[0]);
      (*z) = createCSRMatrix((Sparsity*)z_csr);
      csrMat * cz = getCSRMat(*z);
      for(int xrw = 0; xrw < x_nr; ++xrw)
        for(int ycl = 0; ycl < y_nc; ++ycl)
          for(int inr = 0; inr < x_nc; ++inr)
            (*cz)(xrw,ycl) += (*cx)(xrw,inr) * (*cy)(inr,ycl);
    }
  };
  MatVecMult * getSparseMatVecMult()
  {
    static sparseMatVec * mvm = nullptr;
    if(mvm == nullptr)
      mvm = new sparseMatVec;
    return mvm;
  }
  MatMatMult * getSparseMatMatMult()
  {
    static sparseMatMat * mmm = nullptr;
    if(mmm == nullptr)
      mmm = new sparseMatMat;
    return mmm;
  }

}
