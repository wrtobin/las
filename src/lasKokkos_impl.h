#ifndef LAS_KOKKOS_IMPL_H_
#define LAS_KOKKOS_IMPL_H_
#include <Kokkos_Core.hpp>
#include <Kokkos_StaticCrsGraph.hpp>
namespace las
{
  template < typename Val, class Dev >
  struct kokkosMat
  {
    kokkos::StaticCrsGraph<int, Dev, void int> sparsity;
    kokkos::View<Val*,Dev> values;
  }
  void initKokkosLAS(int * argc, char ** argv[], MPI_Comm cm)
  {
    Kokkos::initialize(argc,argv);
  }
  void finalizeKokkosLAS()
  {
    Kokkos::finalize();
  }
  class kokkosMatBuilder : public LasCreateMat
  {
  public:
    
  };
  class kokkosVecBuilder : public LasCreateVec
  {
    
  };
  class kokkos : public LasOps<kokkos>
  {
  public:
    void _zero(Mat * m)
    { }
  };
}
#endif
