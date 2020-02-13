#include "Kokkos_Core.hpp"
#include "input.hpp"
#include "mpi.hpp"

void applyBCs(struct inputConfig cf, Kokkos::View<double**> &u, class mpiBuffers &m, Kokkos::View<int*> nlft, Kokkos::View<int*> nrht, Kokkos::View<int*> nbot, Kokkos::View<int*> ntop, Kokkos::View<int*> nbak, Kokkos::View<int*> nfrt, Kokkos::View<int*> cell_type, int ncells);

