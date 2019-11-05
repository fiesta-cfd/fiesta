set (CMAKE_CXX_COMPILER "/usr/bin/mpicxx")
set (CMAKE_INSTALL_PREFIX "/home/beromer/src/fiesta-dev/install")

project (Fiesta)
cmake_minimum_required (VERSION 3.15)

set (KOKKOS_ENABLE_CUDA "on")
set (KOKKOS_ENABLE_CUDA_LAMBDA "on")
set (KOKKOS_CUDA_DIR "/opt/cuda")
set (KOKKOS_PATH "/home/beromer/src/fiesta-dev/fiesta/kokkos)

add_subdirectory (kokkos)
#include_directories ("/home/beromer/src/fiesta-dev/kokkos-install/include")
file(GLOB SOURCES "src/*.cpp")


add_executable (fiesta ${SOURCES})
