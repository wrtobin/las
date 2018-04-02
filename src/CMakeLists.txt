set(las_install las.h
                lasComm.h
                lasCSR.h
                lasCSR_impl.h
                lasDebug.h
                lasScalar.h
                lasSparse.h
                lasSparse_impl.h)
set(las_headers lasCSRBuilder.h) #non-public headers
set(las_sources lasComm.cc
                lasCSR.cc
                lasCSRBuilder.cc)

if(${BUILD_SPARSKIT})
  set(las_install ${las_install} lasSparskitExterns.h lasSparskit.h lasSparskit_impl.h)
  set(las_sources ${las_sources} lasSparskitExterns.cc)
  set(link_libs ${link_libs} sparskit)
endif()

# look for petsc first since it might have some dependencies other libraries use (parmetis, zoltan, hdf5)
#find_package(MPI) # should actually make sure we're compiling with mpi support/compilers
find_package(PETSc)
find_package(Scorec)

# handle configuretion based on found packages
if(CUDA_FOUND)
  message(STATUS "Building with CUDA support.")
  set(las_source ${las_sources} lasCuSparse.h lasCuSparse_impl.h)
  set(link_libs ${link_libs} ${CUDA_LIBRARIES})
  set(include_dirs ${include_dirs} ${CUDA_INCLUDE_DIRS})
  message(STATUS "CUDA include directories: ${CUDA_INCLUDE_DIRS}")
endif()
if(SCOREC_FOUND)
  set(las_install ${las_install} lasCSRCore.h)
  set(las_sources ${las_sources} lasCSRCore.cc)
  set(link_libs ${link_libs} ${SCOREC_LIBRARIES})
  set(include_dirs ${include_dirs} ${SCOREC_INCLUDE_DIRS})
endif()
if(PETSC_FOUND)
  set(las_install ${las_install} lasPETSc.h lasPETSc_impl.h)
  set(link_libs ${link_libs} ${PETSC_LIBRARIES})
  set(include_dirs ${include_dirs} ${PETSC_INCLUDE_DIRS})
endif()

set(las_headers ${las_headers} ${las_install})

add_library(las ${las_sources} ${las_headers})
target_compile_features(las PUBLIC cxx_std_11)
set_target_properties(las PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_include_directories(las PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include/las> ${include_dirs})
target_link_libraries(las PUBLIC ${link_libs})

install(TARGETS las LIBRARY DESTINATION lib ARCHIVE DESTINATION lib INCLUDES DESTINATION include/las)
