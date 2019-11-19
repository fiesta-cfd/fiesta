#ifndef RKFUNCTION_H
#define RKFUNCTION_H

#include "input.hpp"
#include "Kokkos_Core.hpp"

class rk_func
{

public:
    rk_func(struct inputConfig &cf_, Kokkos::View<double*> & cd_);

    virtual void compute(const Kokkos::View<double****> & mvar, Kokkos::View<double****> & mdvar) = 0;

protected:
    Kokkos::View<double****> mvar;
    Kokkos::View<double****> mdvar;
    Kokkos::View<double*> mcd;
    struct inputConfig cf;
};

#endif
