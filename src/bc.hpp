#include "Kokkos_Core.hpp"
#include "input.hpp"

void applyBCs(struct inputConfig cf, Kokkos::View<double****> &u);
