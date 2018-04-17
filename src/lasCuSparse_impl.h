#ifndef LAS_CUSPARSE_IMPL_H_
#define LAS_CUSPARSE_IMPL_H_
#include "lasSparse.h"
#include "lasSparse_impl.h"
#include "lasDebug.h"
#include "las.h"
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse.h>
namespace las
{
  // opaque types only used for template specialization
  class cuHost;
  class cuDev;
  typedef cuHost cuDefAlloc;
  typedef csrMat cuMat;
  typedef simpleVec cuVec;
  template <>
  inline void alloc<cuHost>(void ** dat, size_t sz)
  {
    DBG(cudaError_t status =)
      cudaMallocHost(dat,sz);
    assert(status == cudaSuccess && "Allocation of cuda pinned host memory failed.");
  }
  template <>
  inline void dealloc<cuHost>(void ** dat)
  {
    DBG(cudaError_t status =)
      cudaFreeHost(dat);
    assert(status == cudaSuccess && "Deallocation of cuda pinned host memmory failed.");
  }
  template <>
  inline void alloc<cuDev>(void ** dat, size_t sz)
  {
    DBG(cudaError_t status =)
      cudaMalloc(dat,sz);
    assert(status == cudaSuccess && "Allocation of cuda device memory failed.");
  }
  template <>
  inline void dealloc<cuDev>(void ** dat)
  {
    DBG(cudaError_t status =)
      cudaFree(dat);
    assert(status == cudaSuccess && "Deallocation of cuda device memory failed.");
  }
  inline Mat * createCuMat(Sparsity * csr)
  {
    return createCSRMatrix(csr);
  }
  inline void destroyCuMat(Mat * m)
  {
    destroyCSRMatrix(m);
  }
  inline Vec * createCuVec(unsigned n)
  {
    return createVector(n);
  }
  inline void destroyCuVec(Vec * v)
  {
    destroyVector(v);
  }
  template <>
  inline LasOps<cusparse> * getLASOps()
  {
    static cusparse * ops = nullptr;
    if (ops == nullptr)
      ops = new cusparse;
    return ops;
  }
  // currently this only works for square matrices
  class cuMatVecMult : public MatVecMult
  {
  public:
    void exec(Mat * a, Vec * x, Vec * y)
    {
      cuMat * A = getCSRMat(a);
      cuVec * X = getSimpleVec(x);
      cuVec * Y = getSimpleVec(y);
      // this *could* be placed in the lasMultiply constructor if we specified the matrix and vector up-front (or at least the csr structure), if we passed the mat and vecs up-front we could also move data to the device as it is set in the relevant CPU data structures instead of all-at-once
      cublasHandle_t  cublas = nullptr;
      cusparseHandle_t cusparse = nullptr;
      cudaStream_t stream = nullptr;
      cusparseMatDescr_t  mat_desc = nullptr;
      // create cublas/cusparse handles, bind a stream to them
      cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
      cublasCreate(&cublas);
      cublasSetStream(cublas, stream);
      cusparseCreate(&cusparse);
      cusparseSetStream(cusparse, stream);
      // create matrix description
      cusparseCreateMatDescr(&mat_desc);
      cusparseSetMatIndexBase(mat_desc, CUSPARSE_INDEX_BASE_ONE);
      cusparseSetMatType(mat_desc, CUSPARSE_MATRIX_TYPE_GENERAL); // also _SYMMETRIX, _HERMITIAN, _TRIANGULAR
      int * dev_csrRows = nullptr;
      int * dev_csrCols = nullptr;
      scalar * dev_csrVals = nullptr;
      scalar * dev_x = nullptr;
      scalar * dev_y = nullptr;
      CSR * csr = A->getCSR();
      int neq = csr->getNumRows();
      int nnz = csr->getNumNonzero();
      alloc<cuDev>((void**)&dev_csrRows, sizeof(int) * (neq + 1));
      alloc<cuDev>((void**)&dev_csrCols, sizeof(int) * (nnz));
      alloc<cuDev>((void**)&dev_csrVals, sizeof(scalar) * (nnz));
      alloc<cuDev>((void**)&dev_x, sizeof(scalar) * (nnz));
      alloc<cuDev>((void**)&dev_y, sizeof(scalar) * (nnz));
      cudaMemcpy(dev_csrRows, csr->getRows(), sizeof(int) * (neq + 1), cudaMemcpyHostToDevice);
      cudaMemcpy(dev_csrCols, csr->getCols(), sizeof(int) * (nnz), cudaMemcpyHostToDevice);
      cudaMemcpy(dev_csrVals, A->getVals(), sizeof(scalar) * (nnz), cudaMemcpyHostToDevice);
      cudaMemcpy(dev_x, &(*X)[0], sizeof(scalar) * (nnz), cudaMemcpyHostToDevice);
      // perform the matrix-vector multiply (this the the general matrix-vector mult. if we have a specific matrix structure (see cusparseSetMatType) there are better blas operations
      const scalar alpha = 1.0;
      const scalar zero = 0.0;
      // http://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrmv_mergepath
      // y = \alpha \times \mathit{op}(\mathbf{A}) \times x + \beta \times y
      // y = a * op(A) * x + b * y;
      cusparseDcsrmv_mp(cusparse,
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        neq,
        neq,
        nnz,
        &alpha, //a
        mat_desc,
        dev_csrVals,
        dev_csrRows,
        dev_csrCols,
        dev_x,
        &zero, //b
        dev_y);
      cudaMemcpy(&(*Y)[0], dev_y, sizeof(scalar) * (nnz), cudaMemcpyDeviceToHost);
      dealloc<cuDev>((void**)&dev_csrRows);
      dealloc<cuDev>((void**)&dev_csrCols);
      dealloc<cuDev>((void**)&dev_csrVals);
      dealloc<cuDev>((void**)&dev_x);
      dealloc<cuDev>((void**)&dev_y);
      cublasDestroy(cublas);
      cusparseDestroy(cusparse);
      cudaStreamDestroy(stream);
      cusparseDestroyMatDescr(mat_desc);
    }
  };
  inline MatVecMult * createCuMatVecMult()
  {
    return new cuMatVecMult;
  }
  class cu_csrMat_csrMatMult : public MatMatMult
  {
  public:
    void exec(Mat * a, Mat * b, Mat ** c)
    {
      cuMat * A = getCSRMat(a); // m x n
      cuMat * B = getCSRMat(b); // n x k
      CSR * a_csr = A->getCSR();
      CSR * b_csr = B->getCSR();
      int m = a_csr->getNumRows();
      int n = a_csr->getNumCols();
      int k = b_csr->getNumCols();
      assert(m == k);
      cublasHandle_t  cublas = nullptr;
      cusparseHandle_t cusparse = nullptr;
      cudaStream_t stream = nullptr;
      cusparseMatDescr_t a_desc = nullptr;
      cusparseMatDescr_t b_desc = nullptr;
      cusparseMatDescr_t c_desc = nullptr;
      // create cublas/cusparse handles, bind a stream to them
      cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
      cublasCreate(&cublas);
      cublasSetStream(cublas, stream);
      cusparseCreate(&cusparse);
      cusparseSetStream(cusparse, stream);
      // create A matrix description and allocate device memory
      cusparseCreateMatDescr(&a_desc);
      cusparseSetMatIndexBase(a_desc, CUSPARSE_INDEX_BASE_ONE);
      cusparseSetMatType(a_desc, CUSPARSE_MATRIX_TYPE_GENERAL); // also _SYMMETRIC, _HERMITIAN, _TRIANGULAR
      int * dev_acsrRows = nullptr;
      int * dev_acsrCols = nullptr;
      scalar * dev_acsrVals = nullptr;
      int annz = a_csr->getNumNonzero();
      alloc<cuDev>((void**)&dev_acsrRows, sizeof(int) * (m + 1));
      alloc<cuDev>((void**)&dev_acsrCols, sizeof(int) * (annz));
      alloc<cuDev>((void**)&dev_acsrVals, sizeof(scalar) * (annz));
      cusparseCreateMatDescr(&b_desc);
      cusparseSetMatIndexBase(b_desc, CUSPARSE_INDEX_BASE_ONE);
      cusparseSetMatType(b_desc, CUSPARSE_MATRIX_TYPE_GENERAL);
      // create B matrix description and allocate device memory
      int * dev_bcsrRows = nullptr;
      int * dev_bcsrCols = nullptr;
      scalar * dev_bcsrVals = nullptr;
      int bnnz = b_csr->getNumNonzero();
      alloc<cuDev>((void**)&dev_bcsrRows, sizeof(int) * (k + 1));
      alloc<cuDev>((void**)&dev_bcsrCols, sizeof(int) * (bnnz));
      alloc<cuDev>((void**)&dev_bcsrVals, sizeof(scalar) * (bnnz));
      // move A and B to the device
      cudaMemcpy(dev_acsrRows, a_csr->getRows(), sizeof(int) * (m + 1), cudaMemcpyHostToDevice);
      cudaMemcpy(dev_acsrCols, a_csr->getCols(), sizeof(int) * (annz), cudaMemcpyHostToDevice);
      cudaMemcpy(dev_acsrVals, A->getVals(), sizeof(scalar) * (annz), cudaMemcpyHostToDevice);
      cudaMemcpy(dev_bcsrRows, b_csr->getRows(), sizeof(int) * (k + 1), cudaMemcpyHostToDevice);
      cudaMemcpy(dev_bcsrCols, b_csr->getCols(), sizeof(int) * (bnnz), cudaMemcpyHostToDevice);
      cudaMemcpy(dev_bcsrVals, B->getVals(), sizeof(scalar) * (bnnz), cudaMemcpyHostToDevice);
      // generate the c structure
      cusparseCreateMatDescr(&c_desc);
      cusparseSetMatIndexBase(c_desc, CUSPARSE_INDEX_BASE_ONE);
      cusparseSetMatType(c_desc, CUSPARSE_MATRIX_TYPE_GENERAL);
      int * dev_ccsrRows = nullptr;
      int * dev_ccsrCols = nullptr;
      scalar * dev_ccsrVals = nullptr;
      int cnnz = 0;
      int * nnzptr = &cnnz;
      int c_base = 0;
      // calculate cnnz
      alloc<cuDev>((void**)&dev_ccsrRows, sizeof(int) * (m + 1));
      cusparseXcsrgemmNnz(cusparse,
                          CUSPARSE_OPERATION_NON_TRANSPOSE,
                          CUSPARSE_OPERATION_NON_TRANSPOSE,
                          m,n,k,
                          a_desc, annz, dev_acsrRows, dev_acsrCols,
                          b_desc, bnnz, dev_bcsrRows, dev_bcsrCols,
                          c_desc, dev_ccsrRows, nnzptr);
      if(nnzptr != NULL)
        cnnz = *nnzptr;
      else
      {
        cudaMemcpy(&cnnz, dev_ccsrRows+m, sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(&c_base, dev_ccsrRows, sizeof(int), cudaMemcpyDeviceToHost);
        cnnz -= c_base;
      }
      alloc<cuDev>((void**)&dev_ccsrCols, sizeof(int) * cnnz);
      alloc<cuDev>((void**)&dev_ccsrVals, sizeof(scalar) * cnnz);
      cusparseDcsrgemm(cusparse,
                       CUSPARSE_OPERATION_NON_TRANSPOSE,
                       CUSPARSE_OPERATION_NON_TRANSPOSE,
                       m,n,k,
                       a_desc, annz, dev_acsrVals, dev_acsrRows, dev_acsrCols,
                       b_desc, bnnz, dev_bcsrVals, dev_bcsrRows, dev_bcsrCols,
                       c_desc, dev_ccsrVals, dev_ccsrRows, dev_ccsrCols);
      int * hst_ccsrRows = nullptr;
      int * hst_ccsrCols = nullptr;
      scalar * hst_ccsrVals = nullptr;
      alloc<cuHost>((void**)&hst_ccsrRows, sizeof(int) * (m + 1));
      alloc<cuHost>((void**)&hst_ccsrCols, sizeof(int) * cnnz);
      alloc<cuHost>((void**)&hst_ccsrVals, sizeof(int) * cnnz);
      cudaMemcpy(&hst_ccsrRows, dev_ccsrRows, sizeof(int) * (m + 1), cudaMemcpyDeviceToHost);
      cudaMemcpy(&hst_ccsrCols, dev_ccsrCols, sizeof(int) * cnnz, cudaMemcpyDeviceToHost);
      cudaMemcpy(&hst_ccsrVals, dev_ccsrVals, sizeof(scalar) * cnnz, cudaMemcpyDeviceToHost);
      CSR * ccsr = new CSR(m,k,cnnz,hst_ccsrRows,hst_ccsrCols);
      *c = createCSRMatrix((Sparsity*)ccsr);
      dealloc<cuDev>((void**)&dev_acsrRows);
      dealloc<cuDev>((void**)&dev_acsrCols);
      dealloc<cuDev>((void**)&dev_acsrVals);
      dealloc<cuDev>((void**)&dev_bcsrRows);
      dealloc<cuDev>((void**)&dev_bcsrCols);
      dealloc<cuDev>((void**)&dev_bcsrVals);
      dealloc<cuDev>((void**)&dev_ccsrRows);
      dealloc<cuDev>((void**)&dev_ccsrCols);
      dealloc<cuDev>((void**)&dev_ccsrVals);
      cusparseDestroyMatDescr(a_desc);
      cusparseDestroyMatDescr(b_desc);
      cusparseDestroyMatDescr(c_desc);
      cudaStreamDestroy(stream);
      cusparseDestroy(cusparse);
      cublasDestroy(cublas);
    }
  };
  inline MatMatMult * createCuCsrMatMatMult()
  {
    return new cu_csrMat_csrMatMult;
  }
  // only works with square matrices, since CSR only supports square matrices
  class cu_csrMat_dMatMult : public MatMatMult
  {
  public:
    void exec(Mat * a, Mat * b, Mat ** c)
    {
      cuMat * A = getCSRMat(a); // m x n
      cuMat * B = getCSRMat(b); // n x k
      CSR * a_csr = A->getCSR();
      CSR * b_csr = B->getCSR();
      int m = a_csr->getNumRows();
      int k = b_csr->getNumCols();
      assert(m == k);
      (void)c;
    }
  };
}
#endif
