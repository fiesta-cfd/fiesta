#include "fiesta.hpp"
#include "input.hpp"
#include "rkfunction.hpp"
#include "Kokkos_Core.hpp"
#include "noise.hpp"
#include "particle.hpp"

class gen2d_func : public rk_func {

public:
    gen2d_func(struct inputConfig &cf_, Kokkos::View<double*> &cd_);

    void compute();
    void preStep();
    void postStep();
    void preSim();
    void postSim();

    FS2D p;          // Pressure
    FS2D T;          // Temperature
    FS3D tvel;
    FS2D rho;      // Total Density
    FS2D fluxx;  // Weno Fluxes in X direction
    FS2D fluxy;  // Weno Fluxes in Y direction
    FS2D_I noise;          // Pressure
    FS1D cd;
    FS4D metrics;
    FS2D J;
//    Kokkos::View<particleStruct*> particles;
    FSP2D particles;

#ifndef NOMPI
    FS4D ls;
    FS4D lr;
    FS4D rs;
    FS4D rr;
    FS4D bs;
    FS4D br;
    FS4D ts;
    FS4D tr;
    
    FS4DH lsH;
    FS4DH lrH;
    FS4DH rsH;
    FS4DH rrH;
    FS4DH bsH;
    FS4DH brH;
    FS4DH tsH;
    FS4DH trH;
#endif
};
