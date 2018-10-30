#ifndef LAS_SPARSKIT_IMPL_H_
#define LAS_SPARSKIT_IMPL_H_
#include "lasSparskitExterns.h"
#include "lasSparse_impl.h"
#include "lasDebug.h"
#include "lasInline.h"
#include "lasCSRBuilder.h"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
namespace las
{
  typedef csrMat skMat;
  typedef lasVec skVec;
  LAS_INLINE skMat * getSparskitMatrix(Mat * m)
  {
    return reinterpret_cast<skMat*>(m);
  }
  LAS_INLINE skVec * getSparskitVector(Vec * v)
  {
    return reinterpret_cast<skVec*>(v);
  }
  LAS_INLINE ScalarMatMult * createSparskitScalarMatMult() 
  {
    return getSparseScalarMatMult();
  }
  LAS_INLINE MatVecMult * getSparskitMatVecMult() 
  {
    return getSparseMatVecMult();
  }
  LAS_INLINE ScalarMatScalarMatAdd * getSparskitScalarMatScalarMatAdd()
  {
    return getSparseScalarMatScalarMatAdd();
  }
  class SparskitLU : public Solve
  {
  protected:
    SparskitBuffers * bfrs;
    double eps;
    friend class SparskitQuickLU;
  public:
    SparskitLU(SparskitBuffers * b, double e) : bfrs(b), eps(e) {}
    SparskitLU(SparskitLU * s, double e) : bfrs(s->bfrs), eps(e) {}
    virtual void solve(Mat * k, Vec * u, Vec * f);
  };
  // only perform the solve, do not decompose the matrix
  class SparskitQuickLU : public SparskitLU
  {
  public:
    SparskitQuickLU(SparskitBuffers * b, double e) : SparskitLU(b,e) {}
    SparskitQuickLU(SparskitLU * lu, double e) : SparskitLU(lu->bfrs,e) {}
    // the matrix k must have a csr format identical to that used previously in a normal SparskitLU solve
    virtual void solve(Mat * k, Vec * u, Vec * f);
  };
  LAS_INLINE Solve * createSparskitLUSolve(SparskitBuffers * b, double eps)
  {
    return new SparskitLU(b,eps);
  }
  LAS_INLINE Solve * createSparskitLUSolve(Solve * slv, double eps)
  {
    return new SparskitLU(reinterpret_cast<SparskitLU*>(slv),eps);
  }
  LAS_INLINE Solve * createSparskitQuickLUSolve(SparskitBuffers * b, double eps)
  {
    return new SparskitQuickLU(b,eps);
  }
  LAS_INLINE Solve * createSparskitQuickLUSolve(Solve * slv, double eps)
  {
    SparskitLU * skt_slv = reinterpret_cast<SparskitLU*>(slv);
    return new SparskitQuickLU(skt_slv,eps);
  }
  /* note for efficiency purposes when the symmetrix matrix market
   * file is written to the standard output, the numver of nonzeros will be missing
   */
  LAS_INLINE void printSparskitMat(std::ostream & o,
                                   Mat * mi,
                                   PrintType tp,
                                   bool symmetric)
  {
    skMat * m = getSparskitMatrix(mi);
    if (tp == PrintType::full)
    {
      int ndofs = m->getCSR()->getNumRows();
      for (int rr = 0; rr < ndofs; ++rr)
      {
        for (int cc = 0; cc < ndofs; ++cc)
        {
          o << (*m)(rr, cc) << ' ';
        }
        o << '\b' << std::endl;
      }
    }
    else if (tp == PrintType::mmarket)
    {
      CSR * csr = m->getCSR();
      double * vals = m->getVals();
      int nnz = csr->getNumNonzero();
      int numRows = csr->getNumRows();
      int numCols = csr->getNumCols();
      int * rows = csr->getRows();
      int * cols = csr->getCols();
      if (symmetric)
      {
        // print the entries on or below the diagonal
        o << "%%MatrixMarket matrix coordinate real symmetric\n";
        o << "%\n";
        int col=0;
        double val=0;
        int row = 1;
        int col_idx = 0;
        int sym_nnz=0;
        //o << numRows << " " << numCols << " " << sym_nnz << "\n";
        o << numRows << " " << numCols << " ";
        std::streampos start = o.tellp();
        o.width(22); o.fill(' '); // give a buffer so we can fill the number of nonzero w/o issue
        o <<"\n"; // this is so that if we print things to stdout it isn't completely broken
        row=1;
        col_idx=0;
        for (int i = 1; i < numRows + 1; ++i)
        {
          int numRowEntries = rows[i] - rows[i - 1];
          for (int j = 0; j < numRowEntries; ++j)
          {
            assert(col_idx < nnz);
            col = cols[col_idx];
            val = 0;
            if (row >= col)
            {
              val = vals[col_idx];
              if(fabs(val) > 1E-15) {
                o << row << " " << col << " " << std::scientific
                  << std::setprecision(std::numeric_limits<double>::digits10 + 1)
                  << val << "\n";
                ++sym_nnz;
              }
            }
            ++col_idx;
          }
          ++row;
        }
        o.clear();
        o.seekp(start);
        o << sym_nnz;
      }
      else
      {
        o << "%%MatrixMarket matrix coordinate real general\n";
        o << "%\n";
        o << numRows << " " << numCols << " " << nnz << "\n";
        int col;
        double val;
        int row = 1;
        int col_idx = 0;
        for (int i = 1; i < numRows + 1; ++i)
        {
          int numRowEntries = rows[i] - rows[i - 1];
          for (int j = 0; j < numRowEntries; ++j)
          {
            assert(col_idx < nnz);
            col = cols[col_idx];
            val = vals[col_idx];
            if(fabs(val) > 1E-15)
              o << row << " " << col << " " << std::scientific
                << std::setprecision(std::numeric_limits<double>::digits10 + 1)
                << val << "\n";
            ++col_idx;
          }
          ++row;
        }
      }
    }
    else
    {
      std::cerr << "Incorrect Matrix Print type. Skipping matrix printing\n";
    }
  }
  LAS_INLINE double getSparskitMatValue(Mat * k, int rr, int cc)
  {
    skMat * m = getSparskitMatrix(k);
    DBG(int ndofs = m->getCSR()->getNumRows());
    assert(rr < ndofs && rr >= 0);
    assert(cc < ndofs && cc >= 0);
    return (*m)(rr,cc);
  }
  LAS_INLINE void setSparskitMatValue(Mat * k, int rr, int cc, double vl)
  {
    skMat * m = getSparskitMatrix(k);
    //DBG(int ndofs = m->getCSR()->getNumRows());
    assert(rr >= 0);
    assert(cc >= 0);
    if(rr<m->getCSR()->getNumRows() && 
       cc < m->getCSR()->getNumCols())
      (*m)(rr,cc) = vl;
  }
  /*
   * Perform a ilut decomposition folowed by lu solve
   * \param \in k  matrix to solve
   * \param \in f force vector
   * \param \out u vector of unknowns
   * \param lfil the maximum number of elements in each row/column of k
   */
  LAS_INLINE void SparskitLU::solve(Mat * k, Vec * u, Vec * f)
  {
    //DBG(bfrs->zero());
    skMat * mat = getSparskitMatrix(k);
    skVec * uv = getSparskitVector(u);
    skVec * fv = getSparskitVector(f);
    CSR * csr = mat->getCSR();
    int ndofs = csr->getNumRows();
    int bfr_lng = 0;
    int ierr = 0;
    do {
    bfr_lng = bfrs->matrixLength();
    ilut_(&ndofs,
          &(*mat)(0,0),
          csr->getCols(),
          csr->getRows(),
          &ndofs,
          &eps,
          bfrs->matrixBuffer(),
          bfrs->colsBuffer(),
          bfrs->rowsBuffer(),
          &bfr_lng,
          bfrs->doubleWorkBuffer(),
          bfrs->intWorkBuffer(),
          &ierr);
      if(ierr == -2 || ierr == -3) {
        // resizesparskit buffers
        std::cerr<<"Initial buffer size too small resizing"<<std::endl;
        bfrs->resizeMatrixBuffer(int(bfrs->matrixLength()*2.0));
      }
    } while(ierr == -2 || ierr == -3);
    if(ierr != 0)
    {
      std::cerr << "ERROR: ilut_ returned error code " << ierr << std::endl;
      switch(ierr) {
        case -1:
          std::cerr<<"Error. input matrix may be wrong (The elimination process has generated a row in L or U whose length is .gt.  n.)\n";
          break;
        case -2:
          std::cerr<<"The matrix L overflows the array al.\n";
          break;
        case -3:
          std::cerr<<"The matrix U overflows the array alu.\n";
          break;
        case -4:
          std::cerr<<"Illegal value for lfil.\n";
          break;
        case -5:
          std::cerr<<"zero row encountered.\n";
          break;
        default:
          std::cerr<<"zero pivot encountered at step number "<<ierr<<".";
          break;

      }
      return;
    }
    lusol_(&ndofs,
           &(*fv)[0],
           &(*uv)[0],
           bfrs->matrixBuffer(),
           bfrs->colsBuffer(),
           bfrs->rowsBuffer());
  }
  LAS_INLINE void SparskitQuickLU::solve(Mat * k, Vec * u, Vec * f)
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
  LAS_INLINE void readSparskitMatLine(std::stringstream & ss,
                                      int & rw,
                                      int & cl,
                                      double & vl)
  {
    if (ss.peek() == ' ') ss.ignore();
    ss >> rw;
    ss >> cl;
    ss >> vl;
  }
  LAS_INLINE Mat * readSparskitMat(std::istream & in, PrintType pt)
  {
    if (pt != PrintType::mmarket)
    {
      std::cerr << "Reading the requested matrix type is not implemented.\n";
      std::abort();
    }
    std::string line;
    std::getline(in, line);
    if (line.find("coordinate") == std::string::npos)
    {
      std::cerr << "Reading the requested matrix market format is not "
                   "supported. The file must be a coordinate type.\n";
      std::abort();
    }
    if (line.find("real") == std::string::npos)
    {
      std::cerr << "Only reading real values is currently supported\n";
      std::abort();
    }
    bool symmetric = line.find("symmetric") != std::string::npos;
    int numRows, numCols, nnz;
    // skip all comment lines until we reach the first data
    while (std::getline(in, line))
    {
      // skip comment lines
      std::stringstream ss(line);
      if (ss.peek() == '%') continue;
      ss >> numRows >> numCols >> nnz;
      break;
    }
    std::streampos start = in.tellg();
    CSRBuilder csrBuilder(numRows, numCols);
    int row, col;
    double val;
    // reset nnz from the actual file so we account for symmetric matrix
    nnz = 0;
    // create sparsity pattern
    while (std::getline(in, line))
    {
      // skip comment lines
      std::stringstream ss(line);
      if (ss.peek() == '%') continue;
      readSparskitMatLine(ss, row, col, val);
      if (row < col && symmetric) {
        std::cerr<<"Symmetric matrix market format is malformed.\n";
        std::cerr<<"Values should only be on the lower left triangle.\nSkipping entry: ";
        std::cerr<<"Row: "<<row<<" Col: "<<col<<" Value: "<<val<<"\n";
        continue;
      }
      csrBuilder.add(row-1, col-1);
      ++nnz;
      if (symmetric)
      {
        // add the terms below the diagonal to the sparsity pattern
        if (row > col)
        {
          csrBuilder.add(col-1, row-1);
          ++nnz;
        }
      }
    }
    Sparsity * csr = csrBuilder.finalize();
    // since we read to EOF we need to reset the EOF flag
    in.clear();
    // go back to the start of the data
    in.seekg(start);
    LasCreateMat * mb = getMatBuilder<sparskit>(0);
    assert(mb);
    Mat * mat = mb->create(nnz, LAS_IGNORE, csr, MPI_COMM_SELF);
    // read in values
    while (std::getline(in, line))
    {
      // skip comment lines
      std::stringstream ss(line);
      if (ss.peek() == '%') continue;
      readSparskitMatLine(ss, row, col, val);
      if (row < col && symmetric) {
        continue;
      }
      setSparskitMatValue(mat, row-1, col-1, val);
      if (symmetric)
      {
        // add the terms below the diagonal to the sparsity pattern
        if (row > col)
        {
          setSparskitMatValue(mat, col-1, row-1, val);
        }
      }
    }
    return mat;
  }
  LAS_INLINE bool isClose(double a, double b, double rtol, double atol)
  {
    return fabs(a - b) <= std::max(rtol * std::max(fabs(a), fabs(b)), atol);
  }
  LAS_INLINE bool sparskitMatClose(Mat * m1i,
                                   Mat * m2i,
                                   double rtol,
                                   double atol)
  {
    skMat * m1 = getSparskitMatrix(m1i);
    skMat * m2 = getSparskitMatrix(m2i);
    CSR * csr1 = m1->getCSR();
    CSR * csr2 = m2->getCSR();
    int nnz1 = csr1->getNumNonzero();
    int nnz2 = csr2->getNumNonzero();
    if (nnz1 != nnz2) return false;
    int numRows1 = csr1->getNumRows();
    int numRows2 = csr2->getNumRows();
    if (numRows1 != numRows2) return false;
    int numCols1 = csr1->getNumCols();
    int numCols2 = csr2->getNumCols();
    if (numCols1 != numCols2) return false;
    int * rows1 = csr1->getRows();
    int * rows2 = csr2->getRows();
    int * cols1 = csr1->getCols();
    int * cols2 = csr2->getCols();
    double * vals1 = m1->getVals();
    double * vals2 = m2->getVals();
    int col_idx = 0;
    int col1, col2;
    double val1 = 0;
    double val2 = 0;
    for (int i = 1; i < numRows1 + 1; ++i)
    {
      if (rows1[i] != rows2[i]) return false;
      int numRowEntries = rows1[i] - rows1[i - 1];
      for (int j = 0; j < numRowEntries; ++j)
      {
        col1 = cols1[col_idx];
        col2 = cols2[col_idx];
        if (col1 != col2) return false;
        val1 = vals1[col_idx];
        val2 = vals2[col_idx];
        if (!isClose(val1, val2, rtol, atol)) return false;
        ++col_idx;
      }
    }
    return true;
  }
}  // namespace las
#endif
