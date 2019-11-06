FIESTA_PATH = /home/beromer/src/fiesta-dev/fiesta
BUILD_DIR = /home/beromer/src/fiesta-dev/build
PREFIX = /home/beromer/.local
FIESTA_SRC = $(FIESTA_PATH)/src


default: fiesta

fiesta: kokkos lua cgns
	export OMPI_CXX=$(BUILD_DIR)/kokkos-build/bin/nvcc_wrapper && make -j -C $(BUILD_DIR)/fiesta-build -f Makefile.fiestabuild

kokkos:
	make -C $(BUILD_DIR)/kokkos-build -f Makefile.kokkosbuild

lua:
	make -C $(BUILD_DIR)/lua-build linux

cgns:
	cd $(BUILD_DIR)/cgns-build && cmake $(FIESTA_PATH)/CGNS -D CGNS_ENABLE_HDF5:BOOL=ON -D HDF5_NEED_MPI:BOOL=ON -D CGNS_ENABLE_PARALLEL:BOOL=ON -D CMAKE_INSTALL_PREFIX:STRING="$(BUILD_DIR)/cgns-build" -D CGNS_ENABLE_SCOPING:BOOL=ON
	make install -j -C $(BUILD_DIR)/cgns-build

clean:
	cd $(BUILD_DIR)/fiesta-dev && rm -f *.o *.cuda *.host *.rocm

