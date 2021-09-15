if (NOT Fiesta_BUILD_HDF5)
    find_package(HDF5 QUIET COMPONENTS C)
    if (HDF5_FOUND)
        if(HDF5_IS_PARALLEL)
            message(STATUS "FIESTA: Found HDF5 with parallel support.")
        else()
            message(FATAL_ERROR "FIESTA: HDF5 was found but does not have parallel support.  Configure with -DFiesta_BUILD_HDF5 or -DFiesta_BUILD_ALL to build HDF5 with parallel support.")
        endif()
    else()
        message(STATUS "FIESTA: An HDF5 installation was not found.")
        set (Fiesta_BUILD_HDF5   ON CACHE BOOL "Build hdf5"   FORCE)
    endif()
endif()

if (Fiesta_BUILD_HDF5)
    message(STATUS "FIESTA: HDF5 will be built with parallel support.")

    FetchContent_Declare(hdf5
        GIT_REPOSITORY https://github.com/HDFGroup/hdf5
        GIT_TAG hdf5-1_12_1
        )

    set(BUILD_SHARED_LIBS OFF CACHE BOOL "")
    set(BUILD_TESTING OFF CACHE BOOL "")
    set(HDF5_BUILD_EXAMPLES OFF CACHE BOOL "")
    set(HDF5_DISABLE_COMPILER_WARNINGS ON CACHE BOOL "")
    set(HDF5_BUILD_CPP_LIB OFF CACHE BOOL "")
    set(HDF5_GENERATE_HEADERS OFF CACHE BOOL "")

    if (Fiesta_NO_MPI)
        set(HDF5_ENABLE_PARALLEL OFF CACHE BOOL "")
    else()
        set(HDF5_ENABLE_PARALLEL ON CACHE BOOL "")
    endif()

    FetchContent_MakeAvailable(hdf5)
endif()

include_directories(${HDF5_C_INCLUDE_DIRS})
