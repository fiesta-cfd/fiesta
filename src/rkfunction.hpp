#ifndef RKFUNCTION_H
#define RKFUNCTION_H

#include "input.hpp"
#include "Kokkos_Core.hpp"

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

    virtual void compute(const Kokkos::View<double****> & mvar, Kokkos::View<double****> & mdvar) = 0;

protected:
    Kokkos::View<double****> mvar;
    Kokkos::View<double****> mdvar;
    Kokkos::View<double*> mcd;
    struct inputConfig cf;
};

#endif