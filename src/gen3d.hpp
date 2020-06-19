#include "fiesta.hpp"
#include "input.hpp"
#include "rkfunction.hpp"
#include "Kokkos_Core.hpp"

class gen3d_func : public rk_func{

public:
    gen3d_func(struct inputConfig &cf, Kokkos::View<double*> & mcd);
    void compute();
    void preStep();
    void postStep();
    void preSim();
    void postSim();

    FS5D metrics;
    FS3D p;
    FS3D T;
    FS4D tvel;
    FS3D rho;
    FS3D qx;
    FS3D qy;
    FS3D qz;
    FS3D fluxx;
    FS3D fluxy;
    FS3D fluxz;
    FS5D stressx;
    FS5D stressy;
    FS5D stressz;
    FS4D gradRho;
    FS4D cFlux;
    FS6D mFlux;
    FS1D cd;
#ifndef NOMPI
    FS5D ls, lr, rs, rr, bs, br, ts, tr, hs, hr, fs, fr;
    FS5DH lsH, lrH, rsH, rrH, bsH, brH, tsH, trH, hsH, hrH, fsH, frH;
#endif

};
