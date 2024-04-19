if (NOT Fiesta_BUILD_KOKKOS)
    find_package(Kokkos QUIET)
    if (Kokkos_FOUND)
        if("${DEVICE}" IN_LIST Kokkos_DEVICES)
            message(STATUS "FIESTA: Found Kokkos with ${DEVICE} support.")
        else()
            message(FATAL_ERROR "FIESTA: Kokkos was found but does not have ${DEVICE} support.  Configure with -DFiesta_BUILD_KOKKOS or -DFiesta_BUILD_ALL to build Kokkos with support for ${DEVICE}")
        endif()
    else()
        message(STATUS "FIESTA: A Kokkos installation was not found.")
        set (Fiesta_BUILD_KOKKOS ON CACHE BOOL "Build kokkos" FORCE)
    endif()
endif()

if (Fiesta_BUILD_KOKKOS)
    message(STATUS "FIESTA: Kokkos will be built with ${DEVICE} support.")

    FetchContent_Declare(kokkos
        GIT_REPOSITORY https://github.com/kokkos/kokkos.git
        GIT_TAG 3.5.00
        )

    set(BUILD_TESTING OFF CACHE BOOL "")
    set(Kokkos_CXX_STANDARD 17 STRING "")
    if (Fiesta_CUDA)
        set(Kokkos_ENABLE_CUDA ON CACHE BOOL "")
        set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "")
        set(Kokkos_ENABLE_LAUNCH_COMPILER ON CACHE BOOL "")
    elseif(Fiesta_HIP)
        set(Kokkos_ENABLE_HIP ON CACHE BOOL "")
        #set(Kokkos_ARCH_VEGA900 ON CACHE BOOL "")
    elseif(Fiesta_OPENMP)
        set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
    elseif()
          set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "")
    endif()

    FetchContent_MakeAvailable(kokkos)
endif()
