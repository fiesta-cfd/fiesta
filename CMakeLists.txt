#  Copyright 2019-2021 The University of New Mexico
#
#  This file is part of FIESTA.
#  
#  FIESTA is free software: you can redistribute it and/or modify it under the
#  terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option) any
#  later version.
#  
#  FIESTA is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
#  details.
#  
#  You should have received a copy of the GNU Lesser General Public License
#  along with FIESTA.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.10)
include(ExternalProject)

# Switch to nvcc compiler wrapper if building for cuda
if (CUDA)
  set(CMAKE_CXX_COMPILER "${CMAKE_SOURCE_DIR}/kokkos/bin/nvcc_wrapper")
elseif(AMDGPU)
  set(CMAKE_CXX_COMPILER "hipcc")
endif()

# set project name and languages
project(Fiesta CXX C Fortran)

set(CMAKE_CXX_STANDARD 17)

# get git branch and hash to set version, build type and date variables
execute_process(COMMAND bash -c "cd ${CMAKE_CURRENT_LIST_DIR} && git describe --tags --dirty=+"
                OUTPUT_VARIABLE GIT_VER)
execute_process(COMMAND bash -c "cd ${CMAKE_CURRENT_LIST_DIR} && git rev-parse --abbrev-ref HEAD"
                OUTPUT_VARIABLE GIT_BRANCH)
string(STRIP ${GIT_VER} GIT_VER)
string(STRIP ${GIT_BRANCH} GIT_BRANCH)
set(FIESTA_VER "${GIT_BRANCH}-${GIT_VER}")
set_property(GLOBAL APPEND
    PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/.git/index")
if (LITE)
    set(FIESTA_OPTS "LITE")
else()
    set(FIESTA_OPTS "FULL")
endif()
string(TIMESTAMP BUILD_TIME "%Y-%m-%dT%H:%MZ" UTC)

# compile flags
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -pedantic -Wall -Wextra -fdiagnostics-color=always")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")

message("######## ${CMAKE_BUILD_TYPE} ########")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(Kokkos_ENABLE_DEBUG ON CACHE BOOL "")
    set(Kokkos_ENABLE_DEBUG_BOUNDS_CHECK ON CACHE BOOL "")
endif()

# set kokkos options depending on build type
if (CUDA)
  # detect cuda architecture
  find_package(CUDA)
  CUDA_SELECT_NVCC_ARCH_FLAGS(TEST_ARCH Auto)

  # set kokkos architecture flags depending on detected cuda architecture
  if ("${TEST_ARCH_readable}" STREQUAL "sm_35")
    set(Kokkos_ARCH_KEPLER35 ON CACHE BOOL "")
  elseif ("${TEST_ARCH_readable}" STREQUAL "sm_60")
    set(Kokkos_ARCH_PASCAL60 ON CACHE BOOL "")
  elseif ("${TEST_ARCH_readable}" STREQUAL "sm_61")
    set(Kokkos_ARCH_PASCAL61 ON CACHE BOOL "")
  elseif ("${TEST_ARCH_readable}" STREQUAL "sm_70")
    message("TURNING ON VOLTA70")
    set(Kokkos_ARCH_VOLTA70 ON CACHE BOOL "")
  endif()

  # enable other kokkos cuda options and append build type string
  set(Kokkos_ENABLE_CUDA ON CACHE BOOL "")
  set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "")
  set(Kokkos_ENABLE_CXX17 ON CACHE BOOL "")
  set(FIESTA_OPTS "${FIESTA_OPTS}+CUDA")
elseif(AMDGPU)
  set(Kokkos_ENABLE_HIP ON CACHE BOOL "")
  set(Kokkos_ARCH_VEGA900 ON CACHE BOOL "")
  set(FIESTA_OPTS "${FIESTA_OPTS}+AMD")
elseif(OPENMP)
  # set kokkos openmp optins and append build type string
  set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
  set(FIESTA_OPTS "${FIESTA_OPTS}+OPENMP")
endif()

# Configure Kokkos
include_directories(${Kokkos_INCLUDE_DIRS_RET})
add_subdirectory(kokkos)

# Get and Configure Tempi
#find_package(Git QUIET)
#if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
#    option(GIT_SUBMODULE "Check submodules during build" ON)
#    if(GIT_SUBMODULE)
#        message(STATUS "Submodule update")
#        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
#                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
#        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
#            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
#        endif()
#    endif()
#endif()
#
#if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/tempi/CMakeLists.txt")
#    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
#endif()

#if(USE_TEMPI)
#  add_subdirectory(extern/tempi)
#endif(USE_TEMPI)

# Configure Lua
include_directories(lua/src)
add_subdirectory(lua)

## enforce static build
#SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#SET(BUILD_SHARED_LIBS OFF)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")

# If doing a full build, find mpi
if (NOT LITE)
    # MPI
    find_package(MPI REQUIRED)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
endif()

# set install directory for third party components
set(installDir ${CMAKE_CURRENT_BINARY_DIR}/install)

# If doing a full build, find or build hdf5
if (NOT LITE)
    if (NOT BUILD_HDF5)
        find_package(HDF5)
    endif()

    if (HDF5_FOUND AND HDF5_IS_PARALLEL)
        include_directories(SYSTEM ${HDF5_C_INCLUDE_DIRS})
    else()
        ExternalProject_Add(hdf5Parallel
          SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/hdf5
          INSTALL_DIR ${installDir}
          BUILD_BYPRODUCTS ${installDir}/lib/libhdf5.a
                           ${installDir}/lib/libhdf5_hl.a
                           ${installDir}/lib/libhdf5_tools.a
          CMAKE_ARGS -D BUILD_SHARED_LIBS:BOOL=OFF
                     -D BUILD_TESTING:BOOL=OFF
                     -D HDF5_BUILD_EXAMPLES:BOOL=OFF
                     -D HDF5_DISABLE_COMPILER_WARNINGS:BOOL=ON
                     -D HDF5_ENABLE_PARALLEL:BOOL=ON
                     -D HDF5_BUILD_CPP_LIB:BOOL=OFF
                     -D CMAKE_INSTALL_PREFIX:PATH=${installDir}
        )
    endif()
endif()

# include third party headers (might be empty if third party components are provided elsewhere)
include_directories(${installDir}/include)

# set base source files
set(FIESTA_SRCS src/rk.cpp src/fiesta.cpp src/luaReader.cpp src/bc.cpp
     src/status.cpp src/timer.cpp src/input.cpp src/output.cpp
     src/main.cpp src/rkfunction.cpp src/writer.cpp
     src/cart2d.cpp src/cart3d.cpp src/gen2d.cpp src/gen3d.cpp)

add_executable(fiesta ${FIESTA_SRCS})
target_link_libraries(fiesta PRIVATE
                      kokkos
                      lua)

if (LITE)
    # add LITE only build files and link
    target_sources(fiesta PRIVATE src/vtk.cpp)
else(LITE)
    # add additional build files for non-lite builds
    target_sources(fiesta PRIVATE src/mpi.cpp src/hdf.cpp)
    target_link_libraries(fiesta PRIVATE kokkos)
  
    if (HDF5_FOUND AND HDF5_IS_PARALLEL)
        message("### Using System HDF5")
        # link to built hdf5 library
        target_link_libraries(fiesta PRIVATE ${HDF5_C_LIBRARIES})
    else()
        # link to found hdf5 library
        message("### Building HDF5")
        add_dependencies(fiesta hdf5Parallel)
        target_link_libraries(fiesta PRIVATE 
                              ${installDir}/lib/libhdf5.a
                              ${installDir}/lib/libhdf5_hl.a
                              ${installDir}/lib/libhdf5_tools.a
        )
    endif(HDF5_FOUND AND HDF5_IS_PARALLEL)

    if (TEMPI)
        target_link_libraries(fiesta PRIVATE tempi)
    endif (TEMPI)

    target_link_libraries(fiesta PRIVATE ${MPI_CXX_LIBRARIES})
endif(LITE)

# add built strings to source code
if(CUDA)
    target_compile_definitions(fiesta PRIVATE HAVE_CUDA)
elseif(OPENMP)
    target_compile_definitions(fiesta PRIVATE HAVE_OPENMP)
endif()
target_compile_definitions(fiesta PRIVATE FIESTA_VERSION="${FIESTA_VER}")
target_compile_definitions(fiesta PRIVATE FIESTA_OPTIONS="${FIESTA_OPTS}")
target_compile_definitions(fiesta PRIVATE FIESTA_BTIME="${BUILD_TIME}")
if (LITE)
    target_compile_definitions(fiesta PRIVATE NOMPI)
endif()
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(fiesta PRIVATE DEBUG)
endif()

# set install destination
install(TARGETS fiesta
        RUNTIME DESTINATION)

# print build strings
message("Fiesta Version ${FIESTA_VER}")
message("Fiesta Build Options ${FIESTA_OPTS}")
message("Fiesta Build Time ${BUILD_TIME}")
