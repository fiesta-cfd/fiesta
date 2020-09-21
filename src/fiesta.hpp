#ifndef FIESTA2_HPP
#define FIESTA2_HPP

#include "input.hpp"
#include "writer.hpp"
#include "rkfunction.hpp"
#include <iostream>
#ifndef NOMPI
#include "hdf.hpp"
#else
#include "vtk.hpp"
#endif

namespace Fiesta{
    struct inputConfig initialize(int, char **);
    void initializeSimulation(struct inputConfig&, rk_func*);
    void reportTimers(struct inputConfig&, rk_func*);
    void checkIO(struct inputConfig&, rk_func*, int, double);
    void finalize();
};

#endif
