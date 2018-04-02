#include "las_capi.h"
#include "las.h"
typedef las::Mat las_mat;
typedef las::Mat las_vec;
/* LAS_BACKEND = PetscOps, cuOps, skOps */
LasOps<LAS_BACKEND> * las_ops;
LasSolve * las_solver;
LasMultiply * las_multiplier;
void las_init(int * argc, char * argv[], MPI_Comm cm)
{

}
void las_free()
{

}
void las_zero_mat(las_mat * m)
{
  las_ops->zero(m);
}
void las_assemble_mat(las_mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls)
{
  las_ops->assemble(m,cntr,rws,cntn,cls,vls);
}
void las_set_mat(las_mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls)
{
  las_ops->assemble(m,cntr,rws,cntc,cls,vls);
}
void las_zero_vec(las_vec * v)
{
  las_ops->zero(v);
}
void las_assemble_vec(las_vec * v, int cnt, int * rws, scalar * vls)
{
  las_ops->assemble(v,cnt,rws,vls);
}
void las_set_vec(las_vec * v, int cnt, int * rws, scalar * vls)
{
  las_ops->set(v,cnt,rws,vls);
}
scalar las_norm(las_vec * v)
{
  return las_ops->norm(v);
}
scalar las_dot(las_vec * a, las_vec * b)
{
  return las_ops->dot(a,b);
}
void las_axpy(scalar a, las_vec * x, las_vec * y)
{
  las_ops->axpy(a,x,y);
}
void las_get_vec(las_vec * v, scalar ** vls)
{
  las_ops->get(v,*vls);
}
void las_restore_vec(las_vec * v, scalar ** vls)
{
  las_ops->restore(v,vls);
}
void las_solve(las_mat * k, las_vec * u, las_vec * f)
{
  las_solver->solve(k,u,f);
}
void las_multiply(las_mat * a, las_vec * x, las_vec * y)
{
  las_multiplier->exec(a,x,y);
}
