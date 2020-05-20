#include "fiesta.hpp"
#include "input.hpp"
#include "rkfunction.hpp"
#include "Kokkos_Core.hpp"

class hydroc3d_func : public rk_func{

public:
    hydroc3d_func(struct inputConfig &cf, Kokkos::View<double*> & mcd);
    void compute();
    void preStep();
    void postStep();
    void preSim();
    void postSim();

    FS3D p;
    FS3D T;
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

};
