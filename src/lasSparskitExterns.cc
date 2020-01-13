#include "lasSparskitExterns.h"
#include <cassert>
#if defined(BGQ) || defined(__ibmxl__)
void ilut_(int *n,double a[],int ja[],int ia[],int *lfil,double *droptol,double *alu,int *jlu,int *ju,int *iwk,double *w,int *jw,int *ierr)
{
  ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr);
}
void lusol_(int *n,double a[],double *x,double *alu,int *jlu,int *ju)
{
  lusol(n,a,x,alu,jlu,ju);
}
void pgmres_(int *n,int *imk,double *rhs,double sol[],double *vv,double *eps,int *maxits,int *iout,double aa[],int *ja,int *ia,double *alu,int *jlu,int *ju,int *ierr)
{
  pgmres(n,imk,rhs,sol,vv,eps,maxits,iout,aa,ja,ia,alu,jlu,ju,ierr);
}
void amux_(int *n,double x[],double y[],double a[],int ja[],int ia[])
{
  amux(n,x,y,a,ja,ia);
}
#endif
namespace las
{
  SparskitBuffers::SparskitBuffers(int num_dofs)
      : heuristic_length(num_dofs * sqrt(num_dofs) * 100)
      , int_work_array(2 * num_dofs)
      , rows(num_dofs)
      , cols(heuristic_length)
      , double_work_array(num_dofs)
      , matrix(heuristic_length)
  {
  }
  SparskitBuffers::SparskitBuffers(int num_dofs, int num_nonzero)
      : heuristic_length(num_nonzero)
      , int_work_array(2 * num_dofs)
      , rows(num_dofs)
      , cols(heuristic_length)
      , double_work_array(num_dofs)
      , matrix(heuristic_length)
  {
  }
  // note that we use the cols buffer and matrix buffer for items that will be returned
  // in MSR format which means the "correct" heuristic length is NNZ+NDZ+1 where NDZ
  // is the number of diagonal zeros
  SparskitBuffers::SparskitBuffers(int num_dofs, int num_nonzero, int num_diag_zeros)
      : heuristic_length(num_nonzero+num_diag_zeros+1)
      , int_work_array(2 * num_dofs)
      , rows(num_dofs)
      , cols(heuristic_length)
      , double_work_array(num_dofs)
      , matrix(heuristic_length)
  {
  }
  void SparskitBuffers::zero()
  {
    int_work_array.assign(int_work_array.size(),0.0);
    rows.assign(rows.size(),0.0);
    cols.assign(cols.size(),0.0);
    double_work_array.assign(double_work_array.size(),0.0);
    matrix.assign(matrix.size(),0.0);
  }
  void SparskitBuffers::resizeMatrixBuffer(int newSize) {
    // only support increasing the buffer size
    if(newSize > matrixLength())
    {
      heuristic_length = newSize;
      matrix.resize(heuristic_length);
      cols.resize(heuristic_length);
      assert(cols.size() == heuristic_length);
      assert(matrix.size() == heuristic_length);
    }
  }
};
