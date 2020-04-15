#include "fiesta.hpp"
#include "Kokkos_Core.hpp"
#include "input.hpp"
#ifndef NOMPI
#include "mpi.hpp"
#endif

void applyBCs(struct inputConfig cf, FS4D &u, class mpiBuffers &m);
void applyBCs(struct inputConfig cf, FS4D &u);
