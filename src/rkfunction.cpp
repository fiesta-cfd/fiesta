#include "rkfunction.hpp"
#include "Kokkos_Core.hpp"
#include "input.hpp"

rk_func::rk_func(struct inputConfig &cf_, FS1D &cd_)
//rk_func::rk_func(struct inputConfig &cf_, Kokkos::View<double *> &cd_)
    : cf(cf_), mcd(cd_){};
