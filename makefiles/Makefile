FIESTA_PATH = /home/beromer/src/fiesta-dev/fiesta
BUILD_DIR = /home/beromer/src/fiesta-dev/build
PREFIX = /home/beromer/.local
FIESTA_SRC = $(FIESTA_PATH)/src


default: fiesta

fiesta: kokkos
	make -j -C $(BUILD_DIR)/fiesta-build -f Makefile.fiestabuild

kokkos:
	make -C $(BUILD_DIR)/kokkos-build -f Makefile.kokkosbuild

clean:
	cd $(BUILD_DIR)/fiesta-dev && rm -f *.o *.cuda *.host *.rocm

