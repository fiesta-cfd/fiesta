#include "fiesta.hpp"
#include "input.hpp"
#include "rkfunction.hpp"
#include "Kokkos_Core.hpp"
#include "noise.hpp"

class hydro2dvisc_func : public rk_func {

public:
    hydro2dvisc_func(struct inputConfig &cf_, Kokkos::View<double*> &cd_);

    void compute();
    void preStep();
    void postStep();
    void preSim();
    void postSim();

    FS2D p;          // Pressure
    FS2D T;          // Pressure
    FS2D rho;      // Total Density
    FS2D qx;  // Weno Fluxes in X direction
    FS2D qy;  // Weno Fluxes in X direction
    FS2D fluxx;  // Weno Fluxes in X direction
    FS2D fluxy;  // Weno Fluxes in Y direction
    FS4D stressx;  // stress tensor on x faces
    FS4D stressy;  // stress tensor on y faces
    FS3D gradRho;
    FS2D_I noise;          // Pressure
    FS1D cd;
};
