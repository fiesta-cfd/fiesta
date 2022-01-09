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

cmake_minimum_required(VERSION 3.14)
include(ExternalProject)
include(FetchContent)

# set project name and languages
project(Fiesta CXX C Fortran)

# project wide settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default cache variables
set (Fiesta_CUDA         OFF CACHE BOOL "Enable Cuda build"       )
set (Fiesta_HIP          OFF CACHE BOOL "Enable HIP build"        )
set (Fiesta_OPENMP       OFF CACHE BOOL "Enable OpenMP build"     )
set (Fiesta_SERIAL       OFF CACHE BOOL "Enable Serial build"     )
set (Fiesta_BUILD_ALL    OFF CACHE BOOL "Enable super-build"      )
set (Fiesta_BUILD_KOKKOS OFF CACHE BOOL "Build kokkos"            )
set (Fiesta_BUILD_HDF5   OFF CACHE BOOL "Build hdf5"              )
set (Fiesta_BUILD_LUA    OFF CACHE BOOL "Build lua"               )
set (Fiesta_BUILD_FMT    OFF CACHE BOOL "Build fmt"               )
set (Fiesta_BUILD_TESTS  OFF CACHE BOOL "Build all tests"         )
set (Fiesta_ENABLE_DEBUG OFF CACHE BOOL "Enable Fiesta Debugging" )
set (Fiesta_NO_MPI       OFF CACHE BOOL "Build Fiesta without MPI")

# get git branch and hash to set version, build type and date variables
execute_process(COMMAND bash -c "cd ${CMAKE_CURRENT_SOURCE_DIR} && git describe --tags --dirty=+"
                OUTPUT_VARIABLE GIT_VER)
execute_process(COMMAND bash -c "cd ${CMAKE_CURRENT_SOURCE_DIR} && git rev-parse --abbrev-ref HEAD"
                OUTPUT_VARIABLE GIT_BRANCH)
string(STRIP ${GIT_VER} GIT_VER)
string(STRIP ${GIT_BRANCH} GIT_BRANCH)
set(FIESTA_VER "${GIT_BRANCH}-${GIT_VER}")
set_property(GLOBAL APPEND
    PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/.git/index")
string(TIMESTAMP BUILD_TIME "%Y-%m-%dT%H:%MZ" UTC)
if(Fiesta_NO_MPI)
    set(FIESTA_OPTS "NO_MPI")
else()
    set(FIESTA_OPTS "FULL")
endif()

# compile flags
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -pedantic -Wall -Wextra -fdiagnostics-color=always")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")

# set parameters for debug builds
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(Fiesta_ENABLE_DEBUG ON CACHE BOOL "Enable Fiesta Debugging" FORCE)
    set(Kokkos_ENABLE_DEBUG ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_DEBUG_BOUNDS_CHECK ON CACHE BOOL "" FORCE)
endif()

#setup device list
if (Fiesta_CUDA)
    set(DEVICE CUDA)
    set(FIESTA_OPTS "${FIESTA_OPTS}+CUDA")
elseif(Fiesta_HIP)
    set(DEVICE HIP)
    set(FIESTA_OPTS "${FIESTA_OPTS}+HIP")
elseif(Fiesta_OPENMP)
    set(DEVICE OPENMP)
    set(FIESTA_OPTS "${FIESTA_OPTS}+OPENMP")
else()
    set(DEVICE SERIAL)
    set(FIESTA_OPTS "${FIESTA_OPTS}+SERIAL")
    set(Fiesta_SERIAL ON CACHE BOOL "Enable Serial build" FORCE)
endif()

if(Fiesta_BUILD_ALL)
    message(STATUS "FIESTA: Super-build enabled.")
    set (Fiesta_BUILD_KOKKOS ON CACHE BOOL "Build kokkos" FORCE)
    set (Fiesta_BUILD_HDF5   ON CACHE BOOL "Build hdf5"   FORCE)
    set (Fiesta_BUILD_LUA    ON CACHE BOOL "Build lua"    FORCE)
    set (Fiesta_BUILD_FMT    ON CACHE BOOL "Build fmt"    FORCE)
endif()

#find third party components
if (NOT Fiesta_NO_MPI)
    find_package(MPI QUIET REQUIRED)
endif()
include(tpls/kokkos.cmake)
include(tpls/lua.cmake)
include(tpls/fmt.cmake)
include(tpls/hdf5.cmake)
include(tpls/yogrt.cmake)

#import FiestaCore and fiesta targetrs
include_directories(src)
add_subdirectory(src)
add_dependencies(FiestaCore ${YOGRT_PREFIX})

# import tests
if (Fiesta_BUILD_TESTS)
    include(tests/tests.cmake)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set (Fiesta_ENABLE_DEBUG ON CACHE BOOL "Enable Fiesta Debugging" FORCE)
endif()

# print build strings
message("Fiesta Version:       ${FIESTA_VER}" )
message("Fiesta Build Options: ${FIESTA_OPTS}")
message("Fiesta Build Time:    ${BUILD_TIME}" )
