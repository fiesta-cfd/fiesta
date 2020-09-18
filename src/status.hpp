#include "Kokkos_Core.hpp"
#include "fiesta.hpp"
#include "timer.hpp"
#include "rkfunction.hpp"
#include <iomanip>
#include <iostream>

void statusCheck(int cFlag, struct inputConfig cf, rk_func *f, double time,
                                         fiestaTimer &wall, fiestaTimer &sim);
