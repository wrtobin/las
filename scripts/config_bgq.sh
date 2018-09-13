#!/bin/bash -x
PREFIX=/gpfs/u/home/PASC/PASCmrsn/scratch/install/las

  cmake /gpfs/u/home/PASC/PASCmrsn/barn/las/  \
  -DCMAKE_C_COMPILER=`which mpicc` \
  -DCMAKE_CXX_COMPILER=`which mpicxx` \
  -DCMAKE_Fortran_COMPILER=`which mpif77` \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=0 \
  -DSCOREC_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/install/core/lib/cmake/SCOREC/ \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DPETSC_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/install/petsc-3.9.2/xl \
  -DPETSC_ARCH=/ \
  -DBUILD_SPARSKIT=ON \
  -DWITH_KOKKOS=FALSE \
  -DBUILD_TESTS=FALSE
