#include "rkfunction.hpp"
#include "Kokkos_Core.hpp"
#include "input.hpp"

//rk_func::rk_func(struct inputConfig &cf_, FS1D &cd_)
//    : cf(cf_), mcd(cd_){};
rk_func::rk_func(struct inputConfig &cf_) : cf(cf_){};
