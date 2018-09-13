cmake_minimum_required(VERSION 3.0 FATAL_ERROR) # version 3.8 adds cuda as a language, 3.9 has better VS integration
project(las)
set(LAS_MAJOR_VERSION 0)
set(LAS_MINOR_VERSION 1)
set(LAS_PATCH_VERSION 32)
set(LAS_VERSION ${LAS_MAJOR_VERSION}.${LAS_MINOR_VERSION}.${LAS_PATCH_VERSION})
enable_language(CXX)
enable_language(C)
# enable_language(CUDA) re-enable if we need to start compiling .cu kernels
option(BUILD_SPARSKIT "Build sparskit solver library" ON)
if(${BUILD_SPARSKIT})
  enable_language(Fortran)
endif()
cmake_host_system_information(RESULT HOST QUERY HOSTNAME)
message(STATUS "Configuring to build on: ${HOST}")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if("${HOST}" STREQUAL "q.ccni.rpi.edu")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
else()
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so")
endif()

if(${BUILD_SHARED_LIBS})
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
endif()

set(ConfigPackageLocation lib/cmake/las)

#todo : check scalar=double/complex and set flags

# c++11 support now happens through target_compile_features() per-target
# add C++ compiler flags
message(STATUS "Using the ${CMAKE_CXX_COMPILER_ID} compiler")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray")
  # todo :  set cray-specific cxx flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL")
  # todo :  check for regular xl vs bg/q xl compiler
  #         currently assuming bg/q
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qlanglvl=extended0x -O5 -qhot=level=0 -qsimd=auto -qmaxmem=-1 -qstrict -qstrict_induction -qreport")
  # these definitions are BGQ specific, but for now the only place we use XL compiler is on BGQ
  #add_definitions(-DBGQ)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # todo :  set intel-specific cxx flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # todo :  version checks to see which c++ standards are available
  #         and warn/error if the version is incompatible
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # todo :  version checks to see which c++ standards are available
  #         and warn/error if the version is incompatible
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

# add C compiler flags
if("${CMAKE_C_COMPILER_ID}" STREQUAL "Cray")
  # todo :  set cray-specific c flags
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "XL")
  # todo :  set xl-specific c flags
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  # todo :  set intel-specific c flags
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  # todo :  version checks to see which c standards are available
  #         and warn/error if the version is incompatible
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # todo :  version checks to see which c++ standards are available
  #         and warn/error if the version is incompatible
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
endif()

message(STATUS "CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_C_FLAGS}")

# add Fortran compiler flags
if(${BUILD_SPARSKIT})
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Cray")
  # todo : set cray-specific fortran flags
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "XL")
  # todo : set xl-specific fortran flags
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  # todo : set intel-specific fortran flags
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g noarg_temp_created")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -Wall -cpp -fPIC")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -Wall -cpp -fPIC")
endif()

message(STATUS "Configuring sparskit")
add_subdirectory(sparskit)
endif()

option(WITH_KOKKOS "download an build kokkos las kernal" OFF)
if(${WITH_KOKKOS})
  find_package(Git REQUIRED)
  set(kokkos_dir ${CMAKE_CURRENT_SOURCE_DIR}/kokkos)
  if(EXISTS ${kokkos_dir})
    message(STATUS "Existing kokkos directory detected, attempting to `git pull` the newest version.")
    execute_process(COMMAND ${GIT_EXECUTABLE} pull WORKING_DIRECTORY ${kokkos_dir})
  else()
    message(STATUS "No kokkos directory detected, attempting to `git clone` the project.")
    execute_process(COMMAND ${GIT_EXECUTABLE} clone git@github.com:kokkos/kokkos.git ${CMAKE_CURRENT_SOURCE_DIR}/kokkos)
  endif()
  set(KOKKOS_FOUND TRUE)
  add_subdirectory(kokkos)
endif()

message(STATUS "Configuring las")
add_subdirectory(src)

find_package(PETSc) # may influence scorec dependencies (remove eventually)
find_package(SCOREC)
if(SCOREC_FOUND)
  message(STATUS "Configure SCOREC/core integration")
  add_subdirectory(core)
endif()

option(BUILD_TESTS "Build available tests" ON)
if(BUILD_TESTS)
  enable_testing()
  message(STATUS "Configuring tests")
  add_subdirectory(test)
endif(BUILD_TESTS)
