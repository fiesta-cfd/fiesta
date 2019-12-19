#include "input.hpp"
#include "rkfunction.hpp"
#include "Kokkos_Core.hpp"

class wenoc3d_func : public rk_func{

public:
    wenoc3d_func(struct inputConfig &cf, Kokkos::View<double*> & mcd);
    void compute(const Kokkos::View<double****> & mvar, Kokkos::View<int*> nlft, Kokkos::View<int*> nrht, Kokkos::View<int*> nbot, Kokkos::View<int*> ntop, Kokkos::View<int*> nbak, Kokkos::View<int*> nfrt, Kokkos::View<int*> cell_type, Kokkos::View<double****> & mdvar, int ncells);

protected:
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_f;
};
