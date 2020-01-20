#include "Kokkos_Core.hpp"
#include "input.hpp"
#include "mpi.hpp"

void applyBCs(struct inputConfig cf, Kokkos::View<double****> &u, class mpiBuffers &m);
