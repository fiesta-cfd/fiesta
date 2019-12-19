#include "input.hpp"
#include "rkfunction.hpp"
#include "Kokkos_Core.hpp"

class weno2d_func : public rk_func {

public:
    weno2d_func(struct inputConfig &cf_, Kokkos::View<double*> &cd_);

    void compute(const Kokkos::View<double**> & mvar, Kokkos::View<int*> nlft, Kokkos::View<int*> nrht, Kokkos::View<int*> nbot, Kokkos::View<int*> ntop, Kokkos::View<int*> nbak, Kokkos::View<int*> nfrt, Kokkos::View<int*> cell_type, Kokkos::View<double**> & mdvar, int ncells);

private:
    typedef Kokkos::RangePolicy<> policy_f1;
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy_f;
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_f4;
};
