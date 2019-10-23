#ifndef MPI_INIT_H
#define MPI_INIT_H

#include "input.hpp"
#include <mpi.h>
#include "Kokkos_Core.hpp"

struct inputConfig  mpi_init(struct inputConfig cf);
void haloExchange(struct inputConfig cf, Kokkos::View<double****> &deviceV);

#endif
