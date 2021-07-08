if (NOT Fiesta_BUILD_FMT)
    find_package(FMT)
    if(FMT_FOUND)
        message(STATUS "FIESTA: Found FMT.")
    else()
        message(STATUS "FIESTA: A Lua installation was not found.")
        set (Fiesta_BUILD_FMT ON CACHE BOOL "Enable FMT build" FORCE)
    endif()
endif()

if (Fiesta_BUILD_FMT)
    message(STATUS "FIESTA: FMT will be built.")
    FetchContent_Declare(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 7.1.3
     )
    FetchContent_MakeAvailable(fmt)
endif()
