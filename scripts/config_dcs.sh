#!/bin/bash
#mpicc=/gpfs/u/home/PASC/PASCmrsn/scratch/test_compile/mpicc
#mpicxx=/gpfs/u/home/PASC/PASCmrsn/scratch/test_compile/mpicxx
#mpif77=/gpfs/u/home/PASC/PASCmrsn/scratch/test_compile/mpif77
export OMPI_CXX=xlc++_r

  cmake /gpfs/u/home/PASC/PASCmrsn/barn/las/  \
  -DCMAKE_C_COMPILER=`which mpicc` \
  -DCMAKE_CXX_COMPILER=`which mpicxx` \
  -DCMAKE_Fortran_COMPILER=`which mpif77` \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=0 \
  -DSCOREC_DIR=/gpfs/u/home/PASC/PASCmrsn/scratch/dcs/install/core/lib/cmake/SCOREC/ \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/gpfs/u/home/PASC/PASCmrsn/scratch/dcs/install/las  \
  -DPETSC_DIR="$PETSC_DIR" \
  -DPETSC_ARCH="$PETSC_ARCH" \
  -DBUILD_SPARSKIT=ON \
  -DWITH_KOKKOS=FALSE \
  -DBUILD_TESTS=FALSE \
  -DCMAKE_CXX_FLAGS="-Ofast -std=c++11 -Wall" \
  -DCMAKE_C_FLAGS="-Ofast -Wall" \
  -DCMAKE_Fortran_FLAGS="-Ofast -Wall"
  #-DCMAKE_C_FLAGS="-O5" \
  #-DCMAKE_CXX_FLAGS="-O5 -sdt=c++11" \
  #-DCMAKE_Fortran_FLAGS="-O5" \
