KOKKOS_PATH = /users/beromer/src/fiesta-dev/kokkos-install
PREFIX = /users/beromer/opt
FIESTA_SRC = /users/beromer/src/fiesta-dev/fiesta/src
KOKKOS_DEVICES="Cuda"
#KOKKOS_DEVICES="Serial"

KOKKOS_SRC_PATH = ${KOKKOS_PATH}
SRC = $(wildcard $(FIESTA_SRC)/*.cpp)
vpath %.cpp $(sort $(dir $(SRC)))

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = mpicxx
CXXFLAGS = -O3 -llua53 -lm -ldl -lcgns
LINK = ${CXX}
#LDFLAGS = 
EXE = fiesta.cuda
KOKKOS_DEVICES=Cuda
KOKKOS_ARCH=Kepler35
else
CXX = mpicxx
CXXFLAGS = -O3 -llua53 -lm -ldl -lcgns
LINK = ${CXX}
#LDFLAGS =  
EXE = fiesta.serial
KOKKOS_DEVICES=Serial
KOKKOS_ARCH=HSW
endif

OBJ = $(notdir $(SRC:.cpp=.o))
#LIB =

default: build

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE) $(CXXFLAGS)

clean: kokkos-clean 
	rm -f *.o *.cuda *.host *.rocm

install:
	cp fiesta.* $(PREFIX)/bin/.

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(EXTRA_INC) -c $< -o $(notdir $@) $(CXXFLAGS)
