
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
    endif()
endif()

if (Fiesta_BUILD_HDF5)
    message(STATUS "FIESTA: HDF5 will be built with parallel support.")

    set(installDir ${CMAKE_BINARY_DIR}/hdf5Parallel)
    ExternalProject_Add(hdf5Parallel
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/hdf5
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
                   -D HDF5_GENERATE_HEADERS:BOOL=OFF
                   -D CMAKE_INSTALL_PREFIX:PATH=${installDir}
    )
    include_directories(${installDir}/include)
    set(HDF5_C_LIBRARIES "${installDir}/lib/libhdf5.a;${installDir}/lib/libhdf5_hl.a;${installDir}/lib/libhdf5_tools.a")
    add_dependencies(fiestaLib hdf5Parallel)
endif()

include_directories(${HDF5_C_INCLUDE_DIRS})
