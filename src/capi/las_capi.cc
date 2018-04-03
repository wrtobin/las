#include "las_capi.h"
#include "las.h"
#include "lasConfig.h"
/* LASBACKEND = PetscOps, cuOps, skOps */
#define LASBACKEND las::PetscOps
las::LasOps<LASBACKEND> * las_ops;
las::LasSolve * las_solver;
las::LasMultiply * las_multiplier;
las::mat_builder * las_mat_builder;
las::mat_destroyer * las_mat_destroyer;
las::vec_builder * las_vec_builder;
las::vec_destroyer * las_vec_destroyer;
void las_init(int * argc, char * argv[], MPI_Comm cm)
{
  (void)argc;
  (void)argv;
  (void)cm;
}
void las_free()
{

}
las_mat * las_create_mat(unsigned lcl, unsigned gbl, unsigned bs, las_sparsity * sparsity, MPI_Comm cm)
{
  return (las_mat*)(*las_mat_builder)(lcl,gbl,bs,(las::Sparsity*)sparsity,cm);
}
void las_destroy_mat(las_mat * m)
{
  (*las_mat_destroyer)((las::Mat*)m);
}
las_vec * las_create_vec(unsigned lcl, unsigned gbl, unsigned bs, MPI_Comm cm)
{
  return (las_vec*)(*las_vec_builder)(lcl,gbl,bs,cm);
}
void las_destroy_vec(las_vec * v)
{
  (*las_vec_destroyer)((las::Vec*)v);
}
void las_zero_mat(las_mat * m)
{
  las_ops->zero((las::Mat*)m);
}
void las_assemble_mat(las_mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls)
{
  las_ops->assemble((las::Mat*)m,cntr,rws,cntc,cls,vls);
}
void las_set_mat(las_mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls)
{
  las_ops->assemble((las::Mat*)m,cntr,rws,cntc,cls,vls);
}
void las_zero_vec(las_vec * v)
{
  las_ops->zero((las::Vec*)v);
}
void las_assemble_vec(las_vec * v, int cnt, int * rws, scalar * vls)
{
  las_ops->assemble((las::Vec*)v,cnt,rws,vls);
}
void las_set_vec(las_vec * v, int cnt, int * rws, scalar * vls)
{
  las_ops->set((las::Vec*)v,cnt,rws,vls);
}
scalar las_norm(las_vec * v)
{
  return las_ops->norm((las::Vec*)v);
}
scalar las_dot(las_vec * a, las_vec * b)
{
  return las_ops->dot((las::Vec*)a,(las::Vec*)b);
}
void las_axpy(scalar a, las_vec * x, las_vec * y)
{
  las_ops->axpy(a,(las::Vec*)x,(las::Vec*)y);
}
void las_get_vec(las_vec * v, scalar ** vls)
{
  las_ops->get((las::Vec*)v,*vls);
}
void las_restore_vec(las_vec * v, scalar ** vls)
{
  las_ops->restore((las::Vec*)v,*vls);
}
void las_solve(las_mat * k, las_vec * u, las_vec * f)
{
  las_solver->solve((las::Mat*)k,(las::Vec*)u,(las::Vec*)f);
}
void las_multiply(las_mat * a, las_vec * x, las_vec * y)
{
  las_multiplier->exec((las::Mat*)a,(las::Vec*)x,(las::Vec*)y);
}
