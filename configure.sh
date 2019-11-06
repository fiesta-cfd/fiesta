#! /bin/bash

#Defaults
PREFIX=$PWD
DEVICE="Serial"


BUILD_DIR=${PWD}

while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        --fiesta-path*)
        FIESTA_PATH="${key#*=}"
        ;;
        --prefix*)
        PREFIX="${key#*=}"
        ;;
        --device*)
        DEVICE="${key#*=}"
        ;;
        --arch*)
        ARCH="${key#*=}"
        
    esac

    shift
done

if [ "$DEVICE" == "Serial" ]; then
    KOKKOS_CONFIG="--with-serial --arch=${ARCH}"
elif [ "$DEVICE" == "Cuda" ]; then
    KOKKOS_CONFIG="--with-cuda --with-cuda-options=enable_lambda --arch=${ARCH}"
fi

if [ -z "${FIESTA_PATH}" ]; then
    FIESTA_PATH=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
fi

echo -e "BUILD_DIR=${PWD}" >> make.vars
echo -e "FIESTA_PATH=${FIESTA_PATH}" >> make.vars
echo -e "PREFIX=${PREFIX}" >> make.vars
echo -e "KOKKOS_CONFIG=${KOKKOS_CONFIG}" >> make.vars

mkdir -p kokkos-build fiesta-build cgns-build lua-build
cp ${FIESTA_PATH}/makefiles/Makefile.kokkosbuild ${BUILD_DIR}/kokkos-build/.
cp ${FIESTA_PATH}/makefiles/Makefile.fiestabuild ${BUILD_DIR}/fiesta-build/.
cp ${FIESTA_PATH}/makefiles/Makefile ${BUILD_DIR}/.
cp -r ${FIESTA_PATH}/lua/* ${BUILD_DIR}/lua-build/.

