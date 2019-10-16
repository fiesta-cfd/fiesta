KOKKOS_PATH = /users/beromer/Kokkos/kokkos
KOKKOS_SRC_PATH = ${KOKKOS_PATH}
#KOKKOS_DEVICES="OpenMP"
KOKKOS_DEVICES="Cuda"
SRC = $(wildcard /users/beromer/src/loboshok/src/*.cpp)
vpath %.cpp $(sort $(dir $(SRC)))

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
#CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
CXX = mpicxx
CXXFLAGS = -O3 -llua53 -lm -ldl -lcgns --expt-extended-lambda
LINK = ${CXX}
LDFLAGS = 
EXE = unmcode.cuda
KOKKOS_DEVICES = "Cuda,OpenMP"
KOKKOS_ARCH = "SNB,Kepler35"
KOKKOS_CUDA_OPTIONS += "enable_lambda"
else
CXX = mpicxx
#CXX = g++
CXXFLAGS = -O3 -llua53 -lm -ldl -lcgns --expt-extended-lambda
LINK = ${CXX}
LDFLAGS =  
EXE = unmcode.host
KOKKOS_DEVICES = "OpenMP"
KOKKOS_ARCH = "SNB"
endif

DEPFLAGS = -M

OBJ = $(notdir $(SRC:.cpp=.o))
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

test: $(EXE)
	./$(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE) $(CXXFLAGS)

clean: kokkos-clean 
	rm -f *.o *.cuda *.host *.rocm

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(EXTRA_INC) -c $< -o $(notdir $@) $(CXXFLAGS)
	#$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $< -o $(notdir $@)
