cmake_minimum_required(VERSION 3.9 FATAL_ERROR) # version 3.8 adds cuda as a language, 3.9 has better VS integration
project(las)
enable_language(CXX)
enable_language(C)
# enable_language(CUDA) re-enable if we need to start compiling .cu kernels
if(${BUILD_SPARSKIT})
  enable_language(Fortran)
endif()
cmake_host_system_information(RESULT HOST QUERY HOSTNAME)
message(STATUS "Configuring to build on: ${HOST}")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so")
if(BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
endif()

#todo : check scalar=double/complex

# c++11 support now happens through target_compile_features() per-target
# add C++ compiler flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray")
  # todo :  set cray-specific cxx flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL")
  # todo :  check for regular xl vs bg/q xl compiler
  #         currently assuming bg/q
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qflag=w -glanglvl=extended0x")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # todo :  set intel-specific cxx flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
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

find_package(CUDA)

message(STATUS "Configuring las")
add_subdirectory(src)
message(STATUS "Configuring tests")
add_subdirectory(test)
