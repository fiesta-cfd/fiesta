cmake_minimum_required(VERSION 3.10)
include(ExternalProject)

if (CUDA)
  set(CMAKE_CXX_COMPILER "${CMAKE_SOURCE_DIR}/kokkos/bin/nvcc_wrapper")
endif()

project(Fiesta CXX)

execute_process(COMMAND bash -c "cd ${CMAKE_CURRENT_LIST_DIR} && git describe --tags --dirty=+"
                OUTPUT_VARIABLE GIT_VER)
execute_process(COMMAND bash -c "cd ${CMAKE_CURRENT_LIST_DIR} && git rev-parse --abbrev-ref HEAD"
                OUTPUT_VARIABLE GIT_BRANCH)
string(STRIP ${GIT_VER} GIT_VER)
string(STRIP ${GIT_BRANCH} GIT_BRANCH)
set(FIESTA_VER "${GIT_BRANCH}-${GIT_VER}")
if (LITE)
    set(FIESTA_OPTS "LITE")
else()
    set(FIESTA_OPTS "FULL")
endif()
string(TIMESTAMP BUILD_TIME "%Y-%m-%dT%H:%MZ" UTC)

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
  set(FIESTA_OPTS "${FIESTA_OPTS}+CUDA")
elseif(OPENMP)
  set(KOKKOS_ENABLE_OPENMP ON)
  set(FIESTA_OPTS "${FIESTA_OPTS}+OPENMP")
endif()
add_subdirectory(kokkos)

# Kokkos
include_directories(${Kokkos_INCLUDE_DIRS_RET})


# lua
include_directories(lua/src)
add_subdirectory(lua)
if (NOT LITE)
    # MPI
    find_package(MPI REQUIRED)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})

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
endif()

include_directories(${installDir}/include)

# fiesta
set(FIESTA_SRCS src/bc.cpp src/status.cpp src/timer.cpp src/input.cpp src/output.cpp src/main.cpp src/rkfunction.cpp src/cart2d.cpp src/cart3d.cpp)
if (LITE)
    message("Fiesta Sources" ${FIESTA_SRCS})
    add_executable(fiesta ${FIESTA_SRCS})
    target_link_libraries(fiesta
                          kokkos
                          lua
    )
else()
    message("Fiesta Sources" ${FIESTA_SRCS})
    add_executable(fiesta ${FIESTA_SRCS} src/cgns.cpp src/mpi.cpp)
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
endif()

target_compile_definitions(fiesta PRIVATE FIESTA_VERSION="${FIESTA_VER}")
target_compile_definitions(fiesta PRIVATE FIESTA_OPTIONS="${FIESTA_OPTS}")
target_compile_definitions(fiesta PRIVATE FIESTA_BTIME="${BUILD_TIME}")
if (LITE)
    target_compile_definitions(fiesta PRIVATE NOMPI)
endif()


install(TARGETS fiesta
        RUNTIME DESTINATION)

message("Fiesta Version ${FIESTA_VER}")
message("Fiesta Build Options ${FIESTA_OPTS}")
message("Fiesta Build Time ${BUILD_TIME}")
