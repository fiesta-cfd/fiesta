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

set(CMAKE_CXX_STANDARD 17)

# Default cache variables
set (Fiesta_CUDA         OFF CACHE BOOL "Enable Cuda build")
set (Fiesta_HIP          OFF CACHE BOOL "Enable HIP build")
set (Fiesta_OPENMP       OFF CACHE BOOL "Enable OpenMP build")
set (Fiesta_SERIAL       OFF CACHE BOOL "Enable Serial build")
set (Fiesta_BUILD_ALL    OFF CACHE BOOL "Enable super-build")
set (Fiesta_BUILD_KOKKOS OFF CACHE BOOL "Build kokkos")
set (Fiesta_BUILD_HDF5   OFF CACHE BOOL "Build hdf5")
set (Fiesta_BUILD_TESTS  OFF CACHE BOOL "Build all tests")

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

add_library(fiestaLib)

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
endif()

if(Fiesta_BUILD_ALL)
    message(STATUS "FIESTA: Super-build enabled.")
    set (Fiesta_BUILD_KOKKOS ON CACHE BOOL "Build kokkos" FORCE)
    set (Fiesta_BUILD_HDF5   ON CACHE BOOL "Build hdf5" FORCE)
endif()

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
include(cmake/kokkos.cmake)

# Configure Lua
add_subdirectory(lua)

add_subdirectory(fmt)


# If doing a full build, find mpi
find_package(MPI QUIET REQUIRED)
#find_package(HDF5 COMPONENTS C)
include(cmake/hdf5.cmake)

# set install directory for third party components
#set(installDir ${CMAKE_CURRENT_BINARY_DIR}/install)
include_directories(src)

# include third party headers (might be empty if third party components are provided elsewhere)
#include_directories(${installDir}/include)

#set base source files
set(FIESTA_SRCS src/rk.cpp src/fiesta.cpp src/luaReader.cpp src/bc.cpp
     src/status.cpp src/timer.cpp src/input.cpp src/output.cpp src/h5.cpp
     src/rkfunction.cpp src/writer.cpp src/mpi.cpp src/hdf.cpp src/xdmf.cpp src/block.cpp)
     

#add_library(fiestaLib ${FIESTA_SRCS})
target_sources(fiestaLib PRIVATE ${FIESTA_SRCS})

add_executable(fiesta src/main.cpp src/cart2d.cpp src/cart3d.cpp src/gen2d.cpp src/gen3d.cpp)

target_link_libraries(fiestaLib
                      Kokkos::kokkos
                      MPI::MPI_CXX
                      ${HDF5_C_LIBRARIES}
                      Lua::lua
                      fmt::fmt
)

target_link_libraries(fiesta
                      fiestaLib
)

# add built strings to source code
if(CUDA)
    target_compile_definitions(fiestaLib PRIVATE HAVE_CUDA)
elseif(OPENMP)
    target_compile_definitions(fiestaLib PRIVATE HAVE_OPENMP)
endif()
target_compile_definitions(fiestaLib PRIVATE FIESTA_VERSION="${FIESTA_VER}")
target_compile_definitions(fiestaLib PRIVATE FIESTA_OPTIONS="${FIESTA_OPTS}")
target_compile_definitions(fiestaLib PRIVATE FIESTA_BTIME="${BUILD_TIME}")
if (LITE)
    target_compile_definitions(fiestaLib PRIVATE NOMPI)
endif()

#add built strings to source code
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
message("MPI Version: ${MPI_CXX_VERSION}")


if (Fiesta_BUILD_TESTS)
    add_subdirectory(catch2)
    add_executable(tests tests/test.cpp src/cart2d.cpp)
    target_link_libraries(tests PRIVATE Catch2::Catch2WithMain Kokkos::kokkos fiestaLib)

    include(CTest)
    include(catch2/extras/Catch.cmake)
    catch_discover_tests(tests)
endif()
