#include "fiesta.hpp"
#include "Kokkos_Core.hpp"
#include "input.hpp"
#include "mpi.hpp"

void applyBCs(struct inputConfig cf, FS4D &u, class mpiBuffers &m);
