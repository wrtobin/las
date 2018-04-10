PREFIX=$DEVROOT/install/las
  cmake .. \
  -DCMAKE_C_COMPILER=`which mpicc` \
  -DCMAKE_CXX_COMPILER=`which mpicxx` \
  -DCMAKE_Fortran_COMPILER=`which mpif90` \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
  -DSCOREC_LIB_DIR=$DEVROOT/install/core/lib \
  -DSCOREC_INCLUDE_DIR=$DEVROOT/install/core/include \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DPETSC_DIR=$PETSC_DIR \
  -DPETSC_ARCH=$PETSC_ARCH \
  -DBUILD_SPARSKIT=ON
  #-DCUDA_CUDART_LIBRARY= #find_cuda has issues with this sometimes
