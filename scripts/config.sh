#!/bin/bash -x
PREFIX=/gpfs/u/scratch/PASC/shared/install/las
#CUDA_INC_PATH=$DEVROOT/install/cuda/9.1/include \
#CUDA_PATH=$DEVROOT/install/cuda/9.1/ \
BUILD=Debug
#BUILD=RelWithDebugInfo
HOSTNAME=`hostname`
if [ "$HOSTNAME" == "q.ccni.rpi.edu" ]; then
  cmake .. \
  -DCMAKE_C_COMPILER=`which mpicc` \
  -DCMAKE_CXX_COMPILER=`which mpicxx` \
  -DCMAKE_Fortran_COMPILER=`which mpif77` \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=0 \
  -DSCOREC_DIR=/gpfs/u/scratch/PASC/shared/install/core/lib/cmake/SCOREC/ \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DPETSC_DIR=/gpfs/u/scratch/PASC/shared/install/petsc-3.9.2/xl \
  -DPETSC_ARCH=/ \
  -DBUILD_SPARSKIT=ON \
  -DWITH_KOKKOS=FALSE \
  -DBUILD_TESTS=FALSE
else
  cmake .. \
  -DCMAKE_BUILD_TYPE=$BUILD \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif77 \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
  -DSCOREC_DIR=$DEVROOT/install/core/$BUILD/lib/cmake/SCOREC/ \
  -DCMAKE_INSTALL_PREFIX=$DEVROOT/install/las/$BUILD  \
  -DPETSC_DIR="$PETSC_DIR" \
  -DPETSC_ARCH="$PETSC_ARCH" \
  -DBUILD_SPARSKIT=ON \
  -DWITH_KOKKOS=FALSE \
  -DBUILD_TESTS=TRUE
fi

#  -DCUDA_TOOLKIT_ROOT_DIR=$DEVROOT/install/cuda/9.1/bin \
#  -DCUDA_CUDART_LIBRARY=$DEVROOT/install/cuda/9.1/lib64/libcudart_static.a

#  -DENABLE_SIMMETRIX=TRUE \
#  -DSCOREC_LIB_DIR=$DEVROOT/install/core/lib \
#  -DSCOREC_INCLUDE_DIR=$DEVROOT/install/core/include \
