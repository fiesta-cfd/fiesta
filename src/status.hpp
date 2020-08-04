#include "Kokkos_Core.hpp"
#include "fiesta.hpp"
#include "timer.hpp"
#include <iomanip>
#include <iostream>

void statusCheck(int cFlag, struct inputConfig cf, FS4D var, int t, double time,
                 fiestaTimer &wall);
