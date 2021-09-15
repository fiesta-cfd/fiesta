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
set(CMAKE_CXX_EXTENSIONS OFF)

# Default cache variables
set (Fiesta_CUDA         OFF CACHE BOOL "Enable Cuda build"      )
set (Fiesta_HIP          OFF CACHE BOOL "Enable HIP build"       )
set (Fiesta_OPENMP       OFF CACHE BOOL "Enable OpenMP build"    )
set (Fiesta_SERIAL       OFF CACHE BOOL "Enable Serial build"    )
set (Fiesta_BUILD_ALL    OFF CACHE BOOL "Enable super-build"     )
set (Fiesta_BUILD_KOKKOS OFF CACHE BOOL "Build kokkos"           )
set (Fiesta_BUILD_HDF5   OFF CACHE BOOL "Build hdf5"             )
set (Fiesta_BUILD_LUA    OFF CACHE BOOL "Build lua"              )
set (Fiesta_BUILD_FMT    OFF CACHE BOOL "Build fmt"              )
set (Fiesta_BUILD_TESTS  OFF CACHE BOOL "Build all tests"        )
set (Fiesta_ENABLE_DEBUG OFF CACHE BOOL "Enable Fiesta Debugging")
set (Fiesta_LITE         OFF CACHE BOOL "Enable Fiesta Debugging")

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
if(Fiesta_LITE)
    set(FIESTA_OPTS "LITE")
else()
    set(FIESTA_OPTS "FULL")
endif()

#Source and Targets
set(FIESTA_LIB_SOURCES
     src/rk.cpp src/fiesta.cpp src/luaReader.cpp src/bc.cpp src/log2.cpp src/reader.cpp
     src/status.cpp src/timer.cpp src/input.cpp src/output.cpp src/h5.cpp
     src/rkfunction.cpp src/writer.cpp src/xdmf.cpp src/block.cpp
)
if(NOT Fiesta_LITE)
    set(FIESTA_LIB_SOURCES ${FIESTA_LIB_SOURCES} src/mpi.cpp)
endif()

set(FIESTA_SOURCES
    src/main.cpp src/cart2d.cpp src/cart3d.cpp src/gen2d.cpp src/gen3d.cpp
)
     
add_library(FiestaCore STATIC ${FIESTA_LIB_SOURCES})
add_executable(fiesta ${FIESTA_SOURCES})

# compile flags
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -pedantic -Wall -Wextra -fdiagnostics-color=always")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")

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
#find third party components
if (NOT Fiesta_LITE)
    find_package(MPI QUIET REQUIRED)
endif()
include(tpls/kokkos.cmake)
include(tpls/hdf5.cmake)
include(tpls/lua.cmake)
include(tpls/fmt.cmake)

# set install directory for third party components
include_directories(src)

#Link libraries
target_link_libraries(FiestaCore
                      Kokkos::kokkos
                      ${HDF5_C_LIBRARIES}
                      ${LUA_LIBRARIES}
                      fmt::fmt
)

if (NOT Fiesta_LITE)
    target_link_libraries(FiestaCore MPI::MPI_CXX)
endif()

target_link_libraries(fiesta FiestaCore)

# add built strings to source code
target_compile_definitions(fiesta     PRIVATE FIESTA_VERSION="${FIESTA_VER}")
target_compile_definitions(FiestaCore PUBLIC  FIESTA_VERSION="${FIESTA_VER}")
target_compile_definitions(fiesta     PRIVATE FIESTA_OPTIONS="${FIESTA_OPTS}")
target_compile_definitions(FiestaCore PUBLIC  FIESTA_OPTIONS="${FIESTA_OPTS}")
target_compile_definitions(fiesta     PRIVATE FIESTA_BTIME="${BUILD_TIME}")
target_compile_definitions(FiestaCore PUBLIC  FIESTA_BTIME="${BUILD_TIME}")

#add build strings to source code
if(Fiesta_CUDA)
    target_compile_definitions(fiesta     PRIVATE HAVE_CUDA)
    target_compile_definitions(FiestaCore PUBLIC  HAVE_CUDA)
endif()
if(Fiesta_HIP)
    target_compile_definitions(fiesta     PRIVATE HAVE_HIP)
    target_compile_definitions(FiestaCore PUBLIC  HAVE_HIP)
endif()
if(Fiesta_OPENMP)
    target_compile_definitions(fiesta     PRIVATE HAVE_OPENMP)
    target_compile_definitions(FiestaCore PUBLIC  HAVE_OPENMP)
endif()
if(Fiesta_SERIAL)
    target_compile_definitions(fiesta     PRIVATE HAVE_SERIAL)
    target_compile_definitions(FiestaCore PUBLIC  HAVE_SERIAL)
endif()
if(Fiesta_LITE)
else()
    target_compile_definitions(fiesta     PRIVATE HAVE_MPI)
    target_compile_definitions(FiestaCore PUBLIC  HAVE_MPI)
endif()


if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set (Fiesta_ENABLE_DEBUG ON CACHE BOOL "Enable Fiesta Debugging" FORCE)
endif()

if(Fiesta_ENABLE_DEBUG)
    target_compile_definitions(fiesta     PRIVATE DEBUG)
    target_compile_definitions(FiestaCore PUBLIC  DEBUG)
endif()

# set install destination
install(TARGETS fiesta RUNTIME DESTINATION)
install(TARGETS FiestaCore RUNTIME DESTINATION)

# print build strings
message("Fiesta Version:       ${FIESTA_VER}"     )
message("Fiesta Build Options: ${FIESTA_OPTS}"    )
message("Fiesta Build Time:    ${BUILD_TIME}"     )

#
if (Fiesta_BUILD_TESTS)
    add_subdirectory(catch2)
    add_executable(tests tests/test.cpp src/cart3d.cpp)
    if (NOT Fiesta_LITE)
        target_compile_definitions(tests     PRIVATE HAVE_MPI)
    endif()
    target_link_libraries(tests PRIVATE Catch2::Catch2WithMain Kokkos::kokkos FiestaCore)

    include(CTest)
    include(catch2/extras/Catch.cmake)
    catch_discover_tests(tests)

    if (NOT Fiesta_LITE)
        add_executable(halotest_unordered tests/halox_unordered.cpp src/cart3d.cpp)
        target_link_libraries(halotest_unordered PRIVATE Kokkos::kokkos FiestaCore)
        add_test(NAME unordered_host_copy_halo_test COMMAND mpirun --oversubscribe -n 8 ./halotest_unordered 1)
        add_test(NAME unordered_gpu_aware_halo_test COMMAND mpirun --oversubscribe -n 8 ./halotest_unordered 2)
        add_test(NAME unordered_gpu_type_halo_test COMMAND mpirun --oversubscribe -n 8 ./halotest_unordered 3)

        add_executable(halotest_ordered tests/halox_ordered.cpp src/cart3d.cpp)
        target_link_libraries(halotest_ordered PRIVATE Kokkos::kokkos FiestaCore)
        add_test(NAME ordered_host_copy_halo_test COMMAND mpirun --oversubscribe -n 8 ./halotest_ordered 1)
        add_test(NAME ordered_gpu_aware_halo_test COMMAND mpirun --oversubscribe -n 8 ./halotest_ordered 2)
    endif()
endif()
