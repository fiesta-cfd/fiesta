cmake_minimum_required(VERSION 3.10)

if (CUDA)
  set(CMAKE_CXX_COMPILER "${CMAKE_SOURCE_DIR}/kokkos/bin/nvcc_wrapper")
endif()

project(Fiesta CXX)

list(APPEND CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -O3)

# CUDA
if (CUDA)
  set(KOKKOS_ENABLE_CUDA ON)
  set(KOKKOS_ENABLE_CUDA_LAMBDA ON)
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

# cgns
set(CGNS_BUILD_SHARED OFF CACHE BOOL "build cgns static")
set(CGNS_ENABLE_64BIT ON CACHE BOOL "build cgns 64bit")
set(CGNS_ENABLE_HDF5 ON CACHE BOOL "enable hdf5")
set(HDF5_NEED_MPI ON CACHE BOOL "hdf5 mpi")
set(CGNS_ENABLE_PARALLEL ON CACHE BOOL "cgns parallel")
set(CGNS_ENABLE_SCOPING ON CACHE BOOL "cgns scoping")

add_subdirectory(cgns)


# fiesta
set(FIESTA_SRCS src/bc.cpp src/cgns.cpp src/input.cpp src/main.cpp src/mpi.cpp src/rkfunction.cpp src/weno2d.cpp src/wenoc3d.cpp)
add_executable(fiesta ${FIESTA_SRCS})
target_link_libraries(fiesta kokkos ${MPI_CXX_LIBRARIES} cgns lua)
