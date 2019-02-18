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
      scalar * c_vals = cm->getVals();
      int * c_cls = csr->getCols();
      int * c_rws = csr->getRows();
      lasVec * sa = getLASVec(a);
      lasVec * sb = getLASVec(b);
      int nr = csr->getNumRows();
      int nc = csr->getNumCols();
      int la = sa->size();
      int lb = sb->size();
      double val;
      assert(nc == la && "Matrix columns and lhs vector length must match");
      assert(nr == lb && "Matrix rows and rhs vector length must match");
      for (int rw = 0; rw < nr; ++rw)
      {
        val = 0;
        for(int idx = c_rws[rw]-1; idx<c_rws[rw+1]-1; ++idx)
        {
          val += c_vals[idx]*(*sa)[c_cls[idx]-1];
        }
        (*sb)[rw] = val;
      }
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
      CSRBuilder bld(x_nr, x_nr);
      for (int xrw = 0; xrw < x_nr; ++xrw)
        for (int ycl = 0; ycl < y_nc; ++ycl)
          for (int inr = 0; inr < x_nc; ++inr)
            if (((*x_csr)(xrw, inr) >= 0 && (*y_csr)(inr, ycl) >= 0))
              bld.add(xrw, ycl);
      Sparsity * z_csr = bld.finalize();
      (*z) = createCSRMatrix(z_csr, true);
      csrMat * cz = getCSRMat(*z);
      for (int xrw = 0; xrw < x_nr; ++xrw)
        for (int ycl = 0; ycl < y_nc; ++ycl)
          for (int inr = 0; inr < x_nc; ++inr)
            (*cz)(xrw, ycl) += (*cx)(xrw, inr) * (*cy)(inr, ycl);
    }
  };
  class sparseScalarMat : public ScalarMatMult
  {
    public:
    void exec(scalar s, Mat * a, Mat ** c)
    {
      csrMat * ca = getCSRMat(a);
      CSR * a_csr = ca->getCSR();
      scalar * vals = ca->getVals();
      int nnz = a_csr->getNumNonzero();
      // if we want to multiply in place
      if (c == NULL)
      {
        for (int i = 0; i < nnz; ++i)
        {
          vals[i] = s * vals[i];
        }
      }
      else
      {
        CSR * c_csr =
            new CSR(a_csr->getNumRows(), a_csr->getNumCols(),
                    a_csr->getNumNonzero(), a_csr->getRows(), a_csr->getCols());
        (*c) = createCSRMatrix(reinterpret_cast<Sparsity *>(c_csr), true);
        scalar * c_vals = getCSRMat(*c)->getVals();
        for (int i = 0; i < nnz; ++i)
        {
          c_vals[i] = s * vals[i];
        }
      }
    }
  };
  template <class ColInstItr, class ValInstItr>
  static void mergeColumns(int * col1_vals,
                           int ncols1,
                           double * vals1,
                           double s1,
                           int * col2_vals,
                           int ncols2,
                           double * vals2,
                           double s2,
                           ColInstItr col_insrt,
                           ValInstItr val_insrt)
  {
    int c1 = 0;
    int c2 = 0;
    while ((c1 < ncols1) || (c2 < ncols2))
    {
      if (c1 >= ncols1)
      {
        *col_insrt = col2_vals[c2];
        *val_insrt = s2 * vals2[c2];
        ++c2;
        ++col_insrt;
        ++val_insrt;
      }
      else if (c2 >= ncols2)
      {
        *col_insrt = col1_vals[c1];
        *val_insrt = s1 * vals1[c1];
        ++c1;
        ++col_insrt;
        ++val_insrt;
      }
      else
      {
        int c1_val = col1_vals[c1];
        int c2_val = col2_vals[c2];
        if (c1_val < c2_val)
        {
          *col_insrt = c1_val;
          *val_insrt = s1 * vals1[c1];
          ++c1;
          ++col_insrt;
          ++val_insrt;
        }
        else if (c2_val < c1_val)
        {
          *col_insrt = c2_val;
          *val_insrt = s2 * vals2[c2];
          ++c2;
          ++col_insrt;
          ++val_insrt;
        }
        else
        {
          *col_insrt = c1_val;
          *val_insrt = s1 * vals1[c1] + s2 * vals2[c2];
          ++c1;
          ++c2;
          ++col_insrt;
          ++val_insrt;
        }
      }
    }
  }
  class sparseScalarMatScalarMatAdd : public MatMatAdd
  {
    public:
    void exec(scalar s1, Mat * a, scalar s2, Mat * b, Mat ** c)
    {
      if (*c != nullptr) destroyCSRMatrix(*c);
      csrMat * ca = getCSRMat(a);
      CSR * a_csr = ca->getCSR();
      csrMat * cb = getCSRMat(b);
      CSR * b_csr = cb->getCSR();
      assert(a_csr->getNumRows() == b_csr->getNumRows() &&
             a_csr->getNumCols() == b_csr->getNumCols());
      std::vector<int> c_rws(a_csr->getNumRows() + 1);
      c_rws[0] = 1;
      std::vector<int> c_cls;
      c_cls.reserve(std::max(a_csr->getNumNonzero(), b_csr->getNumNonzero()));
      std::vector<double> c_vals;
      c_vals.reserve(std::max(a_csr->getNumNonzero(), b_csr->getNumNonzero()));
      int * a_cols;
      int * b_cols;
      double * a_vals;
      double * b_vals;
      int a_ncols;
      int b_ncols;
      for (int rw = 1; rw < a_csr->getNumRows() + 1; ++rw)
      {
        a_cols = &a_csr->getCols()[a_csr->getRows()[rw - 1] - 1];
        a_vals = &ca->getVals()[a_csr->getRows()[rw - 1] - 1];
        a_ncols = a_csr->getRows()[rw] - a_csr->getRows()[rw - 1];
        b_cols = &b_csr->getCols()[b_csr->getRows()[rw - 1] - 1];
        b_vals = &cb->getVals()[b_csr->getRows()[rw - 1] - 1];
        b_ncols = b_csr->getRows()[rw] - b_csr->getRows()[rw - 1];
        // merge the column and value arrays
        mergeColumns(a_cols, a_ncols, a_vals, s1, b_cols, b_ncols, b_vals, s2,
                     std::back_inserter(c_cls), std::back_inserter(c_vals));
        assert(c_cls.size() == c_vals.size());
        c_rws[rw] = c_cls.size() + 1;
      }
      CSR * c_csr = new CSR(a_csr->getNumRows(), a_csr->getNumCols(),
                            c_cls.size(), c_rws, c_cls);
      *c = reinterpret_cast<Mat *>(new csrMat(c_csr, c_vals, true));
    }
  };
  class sparseVecVecAdd : public VecVecAdd
  {
    void exec(scalar s1, Vec * v1, scalar s2, Vec * v2, Vec *& v3)
    {
      lasVec * sv1 = getLASVec(v1);
      lasVec * sv2 = getLASVec(v2);
      assert(sv1->size() == sv2->size());
      lasVec * sv3;
      if (v3)
      {
          destroyVector(v3);
          sv3 = reinterpret_cast<lasVec *>(createVector(sv1->size()));
      }
      else
      {
        sv3 = reinterpret_cast<lasVec *>(createVector(sv1->size()));
      }
      for (int i = 0; i < sv1->size(); ++i)
      {
        (*sv3)[i] = s1 * (*sv1)[i] + s2 * (*sv2)[i];
      }
      v3 = reinterpret_cast<Vec *>(sv3);
    }
  };
  class sparseMatDiagonal : public MatDiagonal
  {
      void exec(scalar s, Mat * m, Vec *& v)
      {
        csrMat * cm = getCSRMat(m);
        CSR * csr = cm->getCSR();
        lasVec * sv;
        assert(csr->getNumCols() == csr->getNumRows());
        if (v)
        {
            destroyVector(v);
            sv = reinterpret_cast<lasVec *>(createVector(csr->getNumRows()));
        }
        else
        {
          sv = reinterpret_cast<lasVec *>(createVector(csr->getNumRows()));
        }
        for(int i=0; i<csr->getNumRows(); ++i) {
          (*sv)[i] = s*(*cm)(i,i);
        }
        v = reinterpret_cast<Vec *>(sv);
      }
  };
  class sparseMatDiagonalInverse : public MatDiagonalInverse
  {
      void exec(scalar s, Mat * m, Vec *& v)
      {
        csrMat * cm = getCSRMat(m);
        CSR * csr = cm->getCSR();
        lasVec * sv;
        assert(csr->getNumCols() == csr->getNumRows());
        if (v)
        {
            destroyVector(v);
            sv = reinterpret_cast<lasVec *>(createVector(csr->getNumRows()));
        }
        else
        {
          sv = reinterpret_cast<lasVec *>(createVector(csr->getNumRows()));
        }
        for(int i=0; i<csr->getNumRows(); ++i) {
          (*sv)[i] = s/(*cm)(i,i);
        }
        v = reinterpret_cast<Vec *>(sv);
      }
  };
  class sparseHadamardProduct : public HadamardProduct
  {
      void exec(Vec * v1, Vec * v2, Vec * v3)
      {
        lasVec * sv1 = getLASVec(v1);
        lasVec * sv2 = getLASVec(v2);
        lasVec * sv3 = getLASVec(v3);
        assert(sv1->size() == sv2->size());
        assert(sv1->size() == sv3->size());
        for(int i=0; i<sv1->size(); ++i)
        {
          (*sv3)[i] = (*sv1)[i]*(*sv2)[i];
        }
      };
  };
  template <>
  MatVecMult * getMatVecMult<sparse>()
  {
    static sparseMatVec * mvm = nullptr;
    if (mvm == nullptr) mvm = new sparseMatVec;
    return mvm;
  }
  template <>
  MatMatMult * getMatMatMult<sparse>()
  {
    static sparseMatMat * mmm = nullptr;
    if (mmm == nullptr) mmm = new sparseMatMat;
    return mmm;
  }
  template <>
  ScalarMatMult * getScalarMatMult<sparse>()
  {
    static sparseScalarMat * smm = nullptr;
    if (smm == nullptr) smm = new sparseScalarMat;
    return smm;
  }
  template <>
  MatMatAdd * getMatMatAdd<sparse>()
  {
    static sparseScalarMatScalarMatAdd * smsma = nullptr;
    if (smsma == nullptr) smsma = new sparseScalarMatScalarMatAdd;
    return smsma;
  }
  template <>
  VecVecAdd * getVecVecAdd<sparse>()
  {
    static sparseVecVecAdd * vva = nullptr;
    if (vva == nullptr) vva = new sparseVecVecAdd;
    return vva;
  }
  template <>
  MatDiagonal * getMatDiagonal<sparse>()
  {
    static sparseMatDiagonal * dia = nullptr;
    if (dia == nullptr) dia = new sparseMatDiagonal;
    return dia;
  }
  template <>
  MatDiagonalInverse * getMatDiagonalInverse<sparse>()
  {
    static sparseMatDiagonalInverse * dia = nullptr;
    if (dia == nullptr) dia = new sparseMatDiagonalInverse;
    return dia;
  }
  template <>
  HadamardProduct * getHadamardProduct<sparse>()
  {
    static sparseHadamardProduct * hp = nullptr;
    if (hp == nullptr) hp = new sparseHadamardProduct;
    return hp;
  }
}  // namespace las
