#include "input.hpp"
#include "rkfunction.hpp"
#include "Kokkos_Core.hpp"

class hydro2d_func : public rk_func {

public:
    hydro2d_func(struct inputConfig &cf_, Kokkos::View<double*> &cd_);

    void compute(const Kokkos::View<double****> & mvar, Kokkos::View<double****> & mdvar);

private:
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy_f;
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_f4;
};
