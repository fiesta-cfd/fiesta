#include "input.hpp"
#include "Kokkos_Core.hpp"

struct weno_func {

    Kokkos::View<double****> mvar;
    Kokkos::View<double****> mdvar;
    Kokkos::View<double*> mcd;

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_f;
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_f4;
    struct inputConfig cf;

    weno_func(struct inputConfig &cf_, const Kokkos::View<double****> & u_,
              Kokkos::View<double****> & k_, Kokkos::View<double*> & cd_);

    void operator()();
};
