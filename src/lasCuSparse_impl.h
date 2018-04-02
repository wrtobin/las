#ifndef LAS_CUSPARSE_IMPL_H_
#define LAS_CUSPARSE_IMPL_H_
#include "lasSparse.h"
#include "lasSparse_impl.h"
#include "lasDebug.h"
#include "las.h"
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cuSparse.h>
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
    DBG(cudaEffor_t status =)
      cudaFreeHost(dat);
    assert(statis == cudaSuccess && "Deallocation of cuda pinned host memmory failed.");
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
    assert(status == cudaSucess && "Deallocation of cuda device memory failed.");
  }
  inline Mat * createCuMat(CSR * csr)
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
  inline LasOps<cuOps> * initCuSparseOps()
  {
    static cuOps * ops = nullptr;
    if (ops == nullptr)
      ops = new cuOps;
    return ops;
  }
  // currently this only works for square matrices
  class cuMultiply : LasMultiply
  {
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
      double * dev_csrVals = nullptr;
      double * dev_x = nullptr;
      double * dev_y = nullptr;
      CSR * csr = A->getCSR();
      int neq = csr->getNumEqs();
      int nnz = csr->getNumNonzero();
      alloc<cuDefAlloc>((void**)&dev_csrRows, sizeof(int) * (neq + 1));
      alloc<cuDefAlloc>((void**)&dev_csrCols, sizeof(int) * (nnz));
      alloc<cuDefAlloc>((void**)&dev_csrVals, sizeof(double) * (nnz));
      alloc<cuDefAlloc>((void**)&dev_x, sizeof(double) * (nnz));
      alloc<cuDefAlloc>((void**)&dev_y, sizeof(double) * (nnz));
      cudaMemcpy(dev_csrRows, csr->getRows(), sizeof(int) * (neq + 1), cudaMemcpyHostToDevice);
      cudaMemcpy(dev_csrCols, csr->getCols(), sizeof(int) * (nnz), cudaMemcpyHostToDevice);
      cudaMemcpy(dev_csrVals, A->getVals(), sizeof(double) * (nnz), cudaMemcpyHostToDevice);
      cudaMemcpy(dev_x, &(*X)[0], sizeof(double) * (nnz), cudaMemcpyHostToDevice);
      // perform the matrix-vector multiply (this the the general matrix-vector mult. if we have a specific matrix structure (see cusparseSetMatType) there are better blas operations
      const double alpha = 1.0;
      const double zero = 0.0;
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
      cudaMemcpy(&(*Y)[0], dev_y, sizeof(double) * (nnz), cudaMemcpyDeviceToHost);
      dealloc<cuDefAlloc>(dev_csrRows);
      dealloc<cuDefAlloc>(dev_csrCols);
      dealloc<cuDefAlloc>(dev_csrVals);
      dealloc<cuDefAlloc>(dev_x);
      dealloc<cuDefAlloc>(dev_y);
      cublasDestroy(cublas);
      cusparseDestroy(cusparse);
      cudaStreamDestroy(stream);
      cusparseDestroyMatDescr(mat_desc);
    }
  };
}
#endif
