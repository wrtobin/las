set(las_install las.h
                lasAlloc.h
                lasComm.h
                lasCSR.h
                lasCSR_impl.h
                lasCSRBuilder.h
                lasCSRBuilder_impl.h
                lasDebug.h
                lasDense.h
                lasDense_impl.h
                lasInline.h
                lasSparse.h
                lasSparseDense.h
                lasSys.h
                lasNNZ.h
                lasSparse_impl.h
                lasVec_impl.h)
set(las_headers ) #non-public headers
set(las_sources las.cc
                lasComm.cc
                lasCSR.cc
                lasCSRBuilder.cc
                lasDense.cc
                lasSparse.cc
                lasSparseDense.cc)

function(to_include header_list include_string)
  foreach(header IN LISTS ${header_list})
    set(string "${string}#include \"${header}\"\n")
  endforeach()
  set(${include_string} "${${include_string}}\n${string}" PARENT_SCOPE)
endfunction(to_include)

if("${BUILD_SPARSKIT}")
  message(STATUS "Building with sparskit support.")
  # define sparskit_headers for lasconfig.h.in
  set(sparskit_api_headers lasSparskit.h)
  set(sparskit_impl_headers lasSparskitExterns.h lasSparskit_impl.h)
  to_include(sparskit_api_headers package_include_headers)
  set(las_install ${las_install} ${sparskit_api_headers} ${sparskit_impl_headers})
  set(las_sources ${las_sources} lasSparskitExterns.cc)
  set(link_libs ${link_libs} sparskit)
endif()

set(USE_MPI 0)
if ("${WITH_MPI}")
    find_package(MPI REQUIRED)
    set(USE_MPI 1)
    include_directories(SYSTEM ${MPI_INCLUDE_DIRS})
endif()
configure_file(lasComm.h.in "${CMAKE_CURRENT_SOURCE_DIR}/lasComm.h" @ONLY)
#find_package(MPI) # should actually make sure we're compiling with mpi support/compilers

if ("${WITH_PETSC}")
  find_package(PETSc REQUIRED)
  message(STATUS "Building with PETSc support.")
  set(petsc_api_headers lasPETSc.h)
  set(petsc_impl_headers lasPETSc_impl.h)
  to_include(petsc_api_headers package_include_headers)
  set(las_install ${las_install} ${petsc_api_headers} ${petsc_impl_headers})
  set(link_libs ${link_libs} ${PETSC_LIBRARIES})
  set(include_dirs ${include_dirs} ${PETSC_INCLUDE_DIRS})
endif()

# handle configuretion based on found packages
find_package(CUDA)
if(CUDA_FOUND)
  message(STATUS "Building with CUDA support.")
  set(cusparse_api_headers lasCuSparse.h)
  set(cusparse_impl_headers lasCuSparse_impl.h)
  to_include(cusparse_api_headers package_include_headers)
  set(las_install ${las_install} ${cusparse_api_headers} ${cusparse_impl_headers})
  set(link_libs ${link_libs} ${CUDA_LIBRARIES})
  set(include_dirs ${include_dirs} ${CUDA_INCLUDE_DIRS})
endif()

if(KOKKOS_FOUND)
  set(kokkos_api_headers lasKokkos.h)
  set(kokkos_impl_headers lasKokkos_impl.h)
  to_include(kokkos_api_headers package_include_headers)
  set(las_install ${las_install} ${kokkos_api_headers} ${kokkos_impl_headers})
  set(link_libs ${link_libs} kokkos)
#  set(include_dirs ${include_dirs} $<BUILD_INTERFACE:${Kokkos_INCLUDE_DIRS_RET}> $<INSTALL_INTERFACE:include/kokkos>)
endif()

configure_file(lasConfig.h.in "${CMAKE_CURRENT_SOURCE_DIR}/lasConfig.h" @ONLY)
set(las_headers ${las_headers} ${las_install} lasConfig.h)

add_library(las STATIC ${las_sources})
if("${HOST}" STREQUAL "q.ccni.rpi.edu")
  target_compile_definitions(las PUBLIC "-DBGQ")
endif()
if(${CMAKE_VERSION} VERSION_GREATER "3.8.2")
  #target_compile_features(las PUBLIC cxx_std_11)
endif()
if(CUDA_FOUND)
  set_target_properties(las PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()
if(MPI_COMPILE_FLAGS)
  set_target_properties(las PROPERTIES COMPILE_FLAGS ${MPI_CXX_COMPILE_DEFINITIONS})
endif()
if(MPI_LINK_FLAGS)
  set_target_properties(las PROPERTIES LINK_FLAGS ${MPI_CXX_LINK_FLAGS})
endif()
set_target_properties(las PROPERTIES PUBLIC_HEADER "${las_headers}")
if("${HOST}" STREQUAL "q.ccni.rpi.edu")
  set_target_properties(las PROPERTIES LINK_SEARCH_START_STATIC 1)
  set_target_properties(las PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()
set(include_dirs $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>  $<INSTALL_INTERFACE:include>  ${include_dirs})
message(STATUS "Include directories for las: ${include_dirs}")
target_include_directories(las PUBLIC ${include_dirs})
target_link_libraries(las PUBLIC ${link_libs} ${MPI_LIBRARIES})
install(TARGETS las EXPORT lasConfig LIBRARY DESTINATION lib ARCHIVE DESTINATION lib COMPONENT kernel PUBLIC_HEADER DESTINATION include)
set(LAS_COMPONENTS las)

function(add_backend backend headers libraries)
  message(STATUS "Adding CAPI backend ${backend}")
  set(BACKEND ${backend})
  to_include(headers backend_header_includes)
  set(capi_header capi/las_capi_${backend}.h)
  set(capi_source capi/las_capi_${backend}.cc)
  configure_file(capi/las_capi.h ${capi_header})
  to_include(capi_header capi_header_include)
  configure_file(capi/las_capi.cc ${capi_source})
  set(capi_lib las_capi_${backend})
  add_library(${capi_lib} STATIC
    ${CMAKE_CURRENT_BINARY_DIR}/capi/las_capi_${backend}.h
    ${CMAKE_CURRENT_BINARY_DIR}/capi/las_capi_${backend}.cc
    ${headers})
  if(${CMAKE_VERSION} VERSION_GREATER "3.8.2")
    #target_compile_features(${capi_lib} PUBLIC cxx_std_11)
  endif()
  target_include_directories(${capi_lib} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>)
  target_link_libraries(${capi_lib} PUBLIC las ${libraries} ${MPI_LIBRARIES})
  if(CUDA_FOUND)
  set_target_properties(${capi_lib} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()
if(MPI_COMPILE_FLAGS)
  set_target_properties(${capi_lib} PROPERTIES COMPILE_FLAGS ${MPI_COMPILE_DEFINITIONS)
endif()
if(MPI_LINK_FLAGS)
  set_target_properties(${capi_lib} PROPERTIES LINK_FLAGS ${MPI_LINK_FLAGS})
endif()
  set_target_properties(${capi_lib} PROPERTIES PUBLIC_HEADER "${CMAKE_CURRENT_BINARY_DIR}/${capi_header};${las_headers}" )
  install(TARGETS ${capi_lib} EXPORT lasConfig LIBRARY DESTINATION lib ARCHIVE DESTINATION lib COMPONENT capi_${backend} PUBLIC_HEADER DESTINATION include)
  set(LAS_COMPONENTS "${LAS_COMPONENTS} ${capi_lib}")
  unset(backend_header_includes)
  unset(BACKEND)
endfunction()

if("${BUILD_SPARSKIT}")
  add_backend(sparskit "${sparskit_api_headers}" sparskit)
endif()
if(CUDA_FOUND)
  add_backend(cusparse "${cusparse_api_headers}" ${CUDA_LIBRARIES})
endif()
if("${WITH_PETSC}")
  target_compile_definitions(las PUBLIC HAVE_PETSC)
  add_backend(petsc "${petsc_api_headers}" ${PETSC_LIBRARIES})
endif()

include(CMakePackageConfigHelpers)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/las/lasConfigVersion.cmake"
  VERSION ${LAS_VERSION} COMPATIBILITY AnyNewerVersion)

#configure_file(${PROJECT_SOURCE_DIR}/cmake/lasConfig.cmake "${CMAKE_CURRENT_BINARY_DIR}/las/lasConfig.cmake" @ONLY)
install(EXPORT lasConfig DESTINATION lib/cmake)
