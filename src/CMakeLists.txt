set(FIESTA_LIB_SOURCES
     rk.cpp fiesta.cpp luaReader.cpp bc.cpp reader.cpp
     status.cpp input.cpp output.cpp h5.cpp
     rkfunction.cpp writer.cpp xdmf.cpp block.cpp
)
if(NOT Fiesta_NO_MPI)
    set(FIESTA_LIB_SOURCES ${FIESTA_LIB_SOURCES} mpi.cpp)
endif()

set(FIESTA_SOURCES
    main.cpp cart2d.cpp cart3d.cpp gen2d.cpp gen3d.cpp
)
     
add_library(FiestaCore STATIC ${FIESTA_LIB_SOURCES})
add_executable(fiesta ${FIESTA_SOURCES})

#Link libraries
target_link_libraries(FiestaCore
    Kokkos::kokkos
    HDF5::HDF5
    ${LUA_LIBRARIES}
    ${YOGRT_LIBS}
    fmt::fmt
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:-lstdc++fs>
    )

if (NOT Fiesta_NO_MPI)
    target_link_libraries(FiestaCore MPI::MPI_CXX)
endif()

target_link_libraries(fiesta FiestaCore)

target_compile_definitions(fiesta     PRIVATE FIESTA_VERSION="${FIESTA_VER}")
target_compile_definitions(FiestaCore PUBLIC  FIESTA_VERSION="${FIESTA_VER}")
target_compile_definitions(fiesta     PRIVATE FIESTA_OPTIONS="${FIESTA_OPTS}")
target_compile_definitions(FiestaCore PUBLIC  FIESTA_OPTIONS="${FIESTA_OPTS}")
target_compile_definitions(fiesta     PRIVATE FIESTA_BTIME="${BUILD_TIME}")
target_compile_definitions(FiestaCore PUBLIC  FIESTA_BTIME="${BUILD_TIME}")

set(DEV_MACRO "HAVE_${DEVICE}")
target_compile_definitions(fiesta     PRIVATE ${DEV_MACRO})
target_compile_definitions(FiestaCore PUBLIC  ${DEV_MACRO})

if(Fiesta_NO_MPI)
    target_compile_definitions(fiesta     PRIVATE HAVE_NOMPI)
    target_compile_definitions(FiestaCore PUBLIC  HAVE_NOMPI)
else()
    target_compile_definitions(fiesta     PRIVATE HAVE_MPI)
    target_compile_definitions(FiestaCore PUBLIC  HAVE_MPI)
endif()

if(Fiesta_ENABLE_DEBUG)
    target_compile_definitions(fiesta     PRIVATE DEBUG)
    target_compile_definitions(FiestaCore PUBLIC  DEBUG)
endif()

if(Fiesta_SINGLE_PRECISION)
    target_compile_definitions(fiesta     PRIVATE HAVE_SINGLE)
    target_compile_definitions(FiestaCore PUBLIC  HAVE_SINGLE)
endif()

# set install destination
install(TARGETS fiesta RUNTIME DESTINATION)
install(TARGETS FiestaCore RUNTIME DESTINATION)


add_executable(fclean 
    fclean.cpp
    )
target_link_libraries(fclean
    Kokkos::kokkos
    ${LUA_LIBRARIES}
    fmt::fmt
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:-lstdc++fs>
    )
install(TARGETS fclean RUNTIME DESTINATION)
