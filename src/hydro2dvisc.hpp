#include "input.hpp"
#include "rkfunction.hpp"
#include "Kokkos_Core.hpp"

class hydro2dvisc_func : public rk_func {

public:
    hydro2dvisc_func(struct inputConfig &cf_, Kokkos::View<double*> &cd_);

    void compute(const FS4D & mvar, FS4D & mdvar);

private:
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy_f;
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_f4;
};
