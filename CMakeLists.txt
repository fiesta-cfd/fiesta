cmake_minimum_required(VERSION 3.10)
include(ExternalProject)

if (CUDA)
  set(CMAKE_CXX_COMPILER "${CMAKE_SOURCE_DIR}/kokkos/bin/nvcc_wrapper")
endif()

project(Fiesta CXX)

list(APPEND CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -O3)

# CUDA
if (CUDA)
  find_package(CUDA)
  CUDA_SELECT_NVCC_ARCH_FLAGS(TEST_ARCH Auto)
  if ("${TEST_ARCH_readable}" STREQUAL "sm_35")
    set(KOKKOS_ARCH "Kepler35" CACHE STRING "")
  elseif ("${TEST_ARCH_readable}" STREQUAL "sm_60")
    set(KOKKOS_ARCH "Pascal60" CACHE STRING "")
  elseif ("${TEST_ARCH_readable}" STREQUAL "sm_61")
    set(KOKKOS_ARCH "Pascal61" CACHE STRING "")
  elseif ("${TEST_ARCH_readable}" STREQUAL "sm_70")
    set(KOKKOS_ARCH "Volta70" CACHE STRING "")
  endif()
  message(STATUS ${KOKKOS_ARCH})
  set(KOKKOS_ENABLE_CUDA ON)
  set(KOKKOS_ENABLE_CUDA_LAMBDA ON)
elseif(OPENMP)
  set(KOKKOS_ENABLE_OPENMP ON)
endif()
add_subdirectory(kokkos)

# Kokkos
include_directories(${Kokkos_INCLUDE_DIRS_RET})

# MPI
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# lua
include_directories(lua/src)
add_subdirectory(lua)

# hdf5
set(installDir ${CMAKE_CURRENT_BINARY_DIR}/install)
ExternalProject_Add(hdf5Parallel
  SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/hdf5
  INSTALL_DIR ${installDir}
  CMAKE_ARGS -D BUILD_SHARED_LIBS:BOOL=OFF
             -D BUILD_TESTING:BOOL=OFF
             -D HDF5_BUILD_EXAMPLES:BOOL=OFF
             -D HDF5_DISABLE_COMPILER_WARNINGS:BOOL=ON
             -D HDF5_ENABLE_PARALLEL:BOOL=ON
             -D HDF5_BUILD_CPP_LIB:BOOL=OFF
             -D CMAKE_INSTALL_PREFIX:PATH=${installDir}
)

# cgns
ExternalProject_Add(cgnsParallel
  SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/cgns
  INSTALL_DIR ${installDir}
  CMAKE_ARGS -D CGNS_BUILD_SHARED:BOOL=OFF
             -D CGNS_USE_SHARED:BOOL=OFF
             -D CGNS_ENABLE_64BIT:BOOL=ON
             -D CGNS_ENABLE_HDF5:BOOL=ON
             -D HDF5_NEED_MPI:BOOL=ON
             -D CGNS_ENABLE_PARALLEL:BOOL=ON
             -D CGNS_ENABLE_SCOPING:BOOL=ON
             -D CMAKE_INSTALL_PREFIX:PATH=${installDir}
             -D CMAKE_PREFIX_PATH:PATH=${installDir}
)

ExternalProject_Add_StepDependencies(cgnsParallel configure hdf5Parallel)
include_directories(${installDir}/include)

# fiesta
set(FIESTA_SRCS src/bc.cpp src/cgns.cpp src/timer.cpp src/input.cpp src/output.cpp src/main.cpp src/mpi.cpp src/rkfunction.cpp src/cart2d.cpp src/cart3d.cpp)
add_executable(fiesta ${FIESTA_SRCS})
add_dependencies(fiesta cgnsParallel)
target_link_libraries(fiesta
                      kokkos
                      ${MPI_CXX_LIBRARIES}
                      ${installDir}/lib/libcgns.a
                      ${installDir}/lib/libhdf5.a
                      ${installDir}/lib/libhdf5_hl.a
                      ${installDir}/lib/libhdf5_tools.a
                      lua
)
install(TARGETS fiesta
        RUNTIME DESTINATION)
