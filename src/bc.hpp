#include "Kokkos_Core.hpp"
#include "fiesta.hpp"
#include "input.hpp"
#ifndef NOMPI
#include "mpi.hpp"
#endif

void applyBCs(struct inputConfig cf, FS4D &u);
