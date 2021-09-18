include(CTest)
enable_testing()

include_directories("tests/common")
set(TEST_NAMES 
    bc_outflow
    bc_reflective
    bc_hydrostatic
    )

if (NOT Fiesta_NO_MPI)
    foreach(test_name IN LISTS TEST_NAMES)
        add_executable(${test_name} tests/${test_name}.cpp tests/common/test.cpp src/cart3d.cpp)
        target_link_libraries(${test_name} PRIVATE Kokkos::kokkos FiestaCore)
        add_test(NAME ${test_name} COMMAND mpirun -n 1 ./tests/${test_name})
        set_target_properties( ${test_name}
            PROPERTIES
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
            )
    endforeach()

    add_executable(halotest_unordered tests/halox_unordered.cpp src/cart3d.cpp)
    target_link_libraries(halotest_unordered PRIVATE Kokkos::kokkos FiestaCore)
    add_test(NAME halox_unordered_host_copy COMMAND mpirun --oversubscribe -n 8 ./tests/halotest_unordered 1)
    add_test(NAME halox_unordered_gpu_aware COMMAND mpirun --oversubscribe -n 8 ./tests/halotest_unordered 2)
    add_test(NAME halox_unordered_gpu_type COMMAND mpirun --oversubscribe -n 8 ./tests/halotest_unordered 3)

    add_executable(halotest_ordered tests/halox_ordered.cpp src/cart3d.cpp)
    target_link_libraries(halotest_ordered PRIVATE Kokkos::kokkos FiestaCore)
    add_test(NAME halox_ordered_host_copy COMMAND mpirun --oversubscribe -n 8 ./tests/halotest_ordered 1)
    add_test(NAME halox_ordered_gpu_aware COMMAND mpirun --oversubscribe -n 8 ./tests/halotest_ordered 2)

    set_target_properties( halotest_ordered halotest_unordered
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
        )
endif()
