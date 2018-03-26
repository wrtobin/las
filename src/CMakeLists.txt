set(las_install las.h lasComm.h lasCSR.h lasCSR_impl.h lasDebug.h lasSparskitExterns.h lasSparskit.h lasSparskit_impl.h)
set(las_sources lasComm.cc lasCSR.cc lasCSRBuilder.cc lasSparskitExterns.cc lasSparskit.cc)

# look for petsc first since it might have some dependencies other libraries use (parmetis, zoltan, hdf5)
find_package(PETSc)
find_package(Scorec REQUIRED)
set(link_libs ${link_libs} ${SCOREC_LIBRARIES})
set(include_dirs ${include_dirs} ${SCOREC_INCLUDE_DIRS})
if(PETSC_FOUND)
  set(las_install ${las_install} lasPETSc.h lasPETSc_impl.h)
  set(las_sources ${las_sources} lasPETSc.cc)
  set(link_libs ${link_libs} ${PETSC_LIBRARIES})
  set(include_dirs ${include_dirs} ${PETSC_INCLUDE_DIRS})
endif()

set(las_headers ${las_headers} ${las_install})

add_library(las ${las_sources} ${las_headers})
target_include_directories(las PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include/las> ${include_dirs})
target_link_libraries(las PUBLIC sparskit ${link_libs})

install(TARGETS las LIBRARY DESTINATION lib ARCHIVE DESTINATION lib INCLUDES DESTINATION include/las)
