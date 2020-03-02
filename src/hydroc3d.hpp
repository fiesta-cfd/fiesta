#include "fiesta.hpp"
#include "input.hpp"
#include "rkfunction.hpp"
#include "Kokkos_Core.hpp"

class hydroc3d_func : public rk_func{

public:
    hydroc3d_func(struct inputConfig &cf, Kokkos::View<double*> & mcd);
    void compute(const FS4D & mvar, FS4D & mdvar);

protected:
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_f;
};
