#ifndef RKFUNCTION_H
#define RKFUNCTION_H

#include "fiesta.hpp"
#include "input.hpp"
#include "Kokkos_Core.hpp"
#include <map>
#include <string>
#include "timer.hpp"

/***
 *
 * Base class for Runge-Kutta function arguments. This is the main feature
 * enabling modularity.  Various schemes or physics can be implemented in
 * classes derived from this class.  Derived classes must implement a 
 * constructor (which initializes the base constructor) and a compute()
 * function which takes a 4d view as the first argument for input and a
 * 4d view as the second argument for output.
 *
 **/

class rk_func
{

public:
    rk_func(struct inputConfig &cf_, Kokkos::View<double*> & cd_);

    virtual void compute() = 0;
    virtual void preStep() = 0;
    virtual void postStep() = 0;
    virtual void preSim() = 0;
    virtual void postSim() = 0;
    //virtual void compute(const FS4D & mvar, FS4D & mdvar) = 0;
    
    FS4D var;
    FS4D dvar;
    FS4D tmp1;
//    FS4D tmp2;
    FS4D grid;

    std::map<std::string, fiestaTimer> timers;
    //std::map<std::string, int> timers;

protected:
    Kokkos::View<double*> mcd;
    struct inputConfig cf;
};

#endif
