#ifndef HYDRO2D_H
#define HYDRO2D_H

#include "input.hpp"
#include "rkfunction.hpp"
#include "Kokkos_Core.hpp"

class hydro2d_func : public rk_func {

public:
    hydro2d_func(struct inputConfig &cf_, Kokkos::View<double*> &cd_);

    //void compute(const FS4D & mvar, FS4D & mdvar);
    void compute();
    FS2D p;
    FS2D rho;
    FS2D wenoy;
    FS2D wenox;
    FS1D cd;

private:
    //policy_f ghost_pol;
    //policy_f face_pol;
    //policy_f cell_pol;
};

#endif
