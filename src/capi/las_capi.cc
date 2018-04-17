@backend_header_includes@
@capi_header_include@
#include "las.h"
las::Solve * las_solver;
las::MatVecMult * las_mv_multiplier;
las::MatMatMult * las_mm_multiplier;
las::LasCreateMat * las_mat_factory;
las::LasCreateVec * las_vec_factory;
las_ops * las(init)(int * argc, char * argv[], MPI_Comm cm)
{
  (void)argc;
  (void)argv;
  (void)cm;
  return (las_ops*)las::getLASOps<las::BACKEND>();
}
las_ops * las(get_ops)()
{
  return (las_ops*)las::getLASOps<las::BACKEND>();
}
void las(free)()
{

}
las_mat * las(create_mat)(unsigned lcl, unsigned bs, las_sparsity * sparsity, MPI_Comm cm)
{
  return (las_mat*)(*las_mat_factory).create(lcl,bs,(las::Sparsity*)sparsity,cm);
}
void las(destroy_mat)(las_mat * m)
{
  (*las_mat_factory).destroy((las::Mat*)m);
}
las_vec * las(create_vec)(unsigned lcl, unsigned bs, MPI_Comm cm)
{
  return (las_vec*)(*las_vec_factory).create(lcl,bs,cm);
}
void las(destroy_vec)(las_vec * v)
{
  (*las_vec_factory).destroy((las::Vec*)v);
}
void las(zero_mat)(las_ops * las, las_mat * m)
{
  ((las::LasOps<las::BACKEND>*)las)->zero((las::Mat*)m);
}
void las(assemble_mat)(las_ops * las, las_mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls)
{
  ((las::LasOps<las::BACKEND>*)las)->assemble((las::Mat*)m,cntr,rws,cntc,cls,vls);
}
void las(set_mat)(las_ops * las, las_mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls)
{
  ((las::LasOps<las::BACKEND>*)las)->assemble((las::Mat*)m,cntr,rws,cntc,cls,vls);
}
void las(get_mat)(las_ops * las, las_mat * m, int cntr, int * rws, int cntc, int * cls, scalar ** vls)
{
  ((las::LasOps<las::BACKEND>*)las)->get((las::Mat*)m,cntr,rws,cntc,cls,vls);
}
void las(zero_row)(las_ops * las, las_mat * m, int rw)
{
  ((las::LasOps<las::BACKEND>*)las)->zero((las::Mat*)m,rw);
}
void las(zero_vec)(las_ops * las, las_vec * v)
{
  ((las::LasOps<las::BACKEND>*)las)->zero((las::Vec*)v);
}
void las(assemble_vec)(las_ops * las, las_vec * v, int cnt, int * rws, scalar * vls)
{
  ((las::LasOps<las::BACKEND>*)las)->assemble((las::Vec*)v,cnt,rws,vls);
}
void las(set_vec)(las_ops * las, las_vec * v, int cnt, int * rws, scalar * vls)
{
  ((las::LasOps<las::BACKEND>*)las)->set((las::Vec*)v,cnt,rws,vls);
}
void las(get_vec)(las_ops * las, las_vec * v, int cnt, int * rws, scalar ** vls)
{
  ((las::LasOps<las::BACKEND>*)las)->get((las::Vec*)v,cnt,rws,vls);
}
scalar las(norm)(las_ops * las, las_vec * v)
{
  return ((las::LasOps<las::BACKEND>*)las)->norm((las::Vec*)v);
}
scalar las(dot)(las_ops * las, las_vec * a, las_vec * b)
{
  return ((las::LasOps<las::BACKEND>*)las)->dot((las::Vec*)a,(las::Vec*)b);
}
void las(axpy)(las_ops * las, scalar a, las_vec * x, las_vec * y)
{
  ((las::LasOps<las::BACKEND>*)las)->axpy(a,(las::Vec*)x,(las::Vec*)y);
}
void las(get_vec_arr)(las_ops * las, las_vec * v, scalar ** vls)
{
  ((las::LasOps<las::BACKEND>*)las)->get((las::Vec*)v,*vls);
}
void las(restore_vec)(las_ops * las, las_vec * v, scalar ** vls)
{
  ((las::LasOps<las::BACKEND>*)las)->restore((las::Vec*)v,*vls);
}
void las(solve)(las_mat * k, las_vec * u, las_vec * f)
{
  las_solver->solve((las::Mat*)k,(las::Vec*)u,(las::Vec*)f);
}
void las(mv_multiply)(las_mat * a, las_vec * x, las_vec * y)
{
  las_mv_multiplier->exec((las::Mat*)a,(las::Vec*)x,(las::Vec*)y);
}
void las(mm_multiply)(las_mat * a, las_mat * b, las_mat ** c)
{
  las_mm_multiplier->exec((las::Mat*)a,(las::Mat*)b,(las::Mat**)c);
}
