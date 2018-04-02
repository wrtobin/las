#ifndef LAS_CAPI_H_
#define LAS_CAPI_H_
#include "lasScalar.h"
#ifdef __cplusplus
extern "C"
{
#endif
/* opaque c types */
struct las_mat;
struct las_vec;
/* utility */
void las_init(int * argc, char * argv[], MPI_Comm cm);
void las_free();
/* mat operations */
void las_zero_mat(las_mat * m);
void las_assemble_mat(las_mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls);
void las_set_mat(las_mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls);
/* vec operations */
void las_zero_vec(las_vec * v);
void las_assemble_vec(las_vec * v, int cnt, int * rws, scalar * vls);
void las_set_vec(las_vec * v, int cnt, int * rws, scalar * vls);
/* arithmetic ops */
scalar las_norm(las_vec * v);
scalar las_dot(las_vec * a, las_vec * b);
void las_axpy(scalar a, las_vec * x, las_vec * y);
/* memory ops */
void las_get_vec(las_vec * v, scalar ** vls);
void las_restore_vec(las_vec * v, scalar ** vls);
/* linear system multiplication and solves */
void las_solve(las_mat * k, las_vec * u, las_vec * f);
void las_mutiply(las_mat * a, las_vec * x, las_vec * y);
#ifdef __cplusplus
}
#endif
#endif
