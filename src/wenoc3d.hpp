#include "input.hpp"
#include "rkfunction.hpp"
#include "Kokkos_Core.hpp"

class wenoc3d_func : public rk_func{

public:
    wenoc3d_func(struct inputConfig &cf, Kokkos::View<double*> & mcd);
    void compute(const Kokkos::View<double****> & mvar, Kokkos::View<double****> & mdvar);

protected:
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_f;
};
