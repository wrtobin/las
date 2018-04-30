#ifndef LAS_CAPI_H_
#define LAS_CAPI_H_
#include "lasSys.h"
#include "lasComm.h"
#define CONCAT_FUNC(pkg,fnc,bck) pkg ## _ ## fnc ## _ ## bck
#define MAKE_PKG_FUNC(pkg,fnc,bck) CONCAT_FUNC(pkg,fnc,bck)
#define las(FNC) MAKE_PKG_FUNC(las,FNC,BACKEND)
#define BACKEND @BACKEND@
#ifdef __cplusplus
extern "C"
{
#endif
/* opaque c types */
struct las_ops;
struct las_mat;
struct las_vec;
struct las_sparsity;
/* utility */
las_ops * las(init)(int * argc, char * argv[], MPI_Comm cm);
void las(free)();
las_ops * las(get_ops)();
/* creation/deletion */
/* las_sparsity * las_create_csr(mesh/field/num) */
/* las_sparsity * las_create_nnz(mesh/field/num) */
las_mat * las(create_mat)(las_ops * las, unsigned lcl, unsigned bs, las_sparsity * sparsity, MPI_Comm cm);
void las(destroy_mat)(las_ops * las, las_mat * m);
las_vec * las(create_vec)(las_ops * las, unsigned sz);
void las(destroy_vec)(las_ops * las, las_vec * v);
/* mat operations */
void las(zero_mat)(las_ops * las, las_mat * m);
void las(assemble_mat)(las_ops * las, las_mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls);
void las(set_mat)(las_ops * las, las_mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls);
void las(get_mat)(las_ops * las, las_mat * m, int cntr, int * rws, int cntc, int * cls, scalar ** vls);
void las(zero_row)(las_ops * las, las_mat * m, int rw);
/* vec operations */
void las(zero_vec)(las_ops * las, las_vec * v);
void las(assemble_vec)(las_ops * las, las_vec * v, int cnt, int * rws, scalar * vls);
void las(set_vec)(las_ops * las, las_vec * v, int cnt, int * rws, scalar * vls);
void las(get_vec)(las_ops * las, las_vec * v, int cnt, int * rws, scalar ** vls);
/* arithmetic ops */
scalar las(norm)(las_ops * las, las_vec * v);
scalar las(dot)(las_ops * las, las_vec * a, las_vec * b);
void las(axpy)(las_ops * las, scalar a, las_vec * x, las_vec * y);
/* memory ops */
void las(get_vec_arr)(las_ops * las, las_vec * v, scalar ** vls);
void las(restore_vec)(las_ops * las, las_vec * v, scalar ** vls);
/* linear system multiplication and solves */
void las(solve)(las_mat * k, las_vec * u, las_vec * f);
void las(multiply)(las_mat * a, las_vec * x, las_vec * y);
#ifdef __cplusplus
}
#endif
#endif
