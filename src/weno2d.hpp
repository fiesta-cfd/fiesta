#include "input.hpp"
#include "Kokkos_Core.hpp"

struct weno2d_func {

    Kokkos::View<double****> mvar;
    Kokkos::View<double****> mdvar;
    Kokkos::View<double*> mcd;

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy_f2;
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_fv2;
    struct inputConfig cf;

    weno2d_func(struct inputConfig &cf_, const Kokkos::View<double****> & u_,
              Kokkos::View<double****> & k_, Kokkos::View<double*> & cd_);

    void operator()();
};
