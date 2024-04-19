set(YOGRT_PREFIX libygrt)
set(YOGRT_URL ${CMAKE_CURRENT_SOURCE_DIR}/tpls/yogrt/libyogrt-1.24.tar)
set(YOGRT_MD5 5cd833467c8746c0eee0294f986daa88)
set(YOGRT_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/${YOGRT_PREFIX})

set(YOGRT_CONFIGURE "--prefix=${YOGRT_INSTALL} --enable-static")
set(YOGRT_INCLUDE ${YOGRT_INSTALL}/include)
set(YOGRT_LIBS ${YOGRT_INSTALL}/lib/yogrt.a)

set(YOGRT_LIBS ${YOGRT_INSTALL}/lib/libyogrt.a)

ExternalProject_Add(${YOGRT_PREFIX}
    PREFIX ${YOGRT_PREFIX}
    URL ${YOGRT_URL}
    URL_MD5 ${YOGRT_MD5}
    CONFIGURE_COMMAND sh -c "./bootstrap && ./configure ${YOGRT_CONFIGURE}"
    BUILD_COMMAND make -j yogrt_build_prefix=${YOGRT_PREFIX}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1

    BUILD_BYPRODUCTS ${YOGRT_LIBS}
)
message(STATUS"YOGRT LIBS: " ${YOGRT_LIBS})
include_directories(${YOGRT_INCLUDE})
