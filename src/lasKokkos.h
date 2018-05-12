#ifndef LAS_KOKKOS_H_
#define LAS_KOKKOS_H_
namespace las
{
  class kokkos;
  template <>
  LasCreateMat * getMatBuilder<kokkos>(int id);
  template <>
  LasCreateMat * getVecBuilder<kokkos>(int id);
  template <>
  LasOps<kokkos> * getLASOps<kokkos>();

}
#endif
