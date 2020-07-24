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

    FS3D p;		// Pressure
    FS3D T;             // Temperature
    FS3D rho;           // Total Density
    FS3D qx;            // Heat Fluxes in X direction
    FS3D qy;            // Heat Fluxes in Y direction  
    FS3D qz;            // Heat Fluxes in Z direction
    FS3D fluxx;         // Weno Fluxes in X direction
    FS3D fluxy;         // Weno Fluxes in Y direction
    FS3D fluxz;         // Weno Fluxes in Z direction
    FS5D stressx;       // Stress tensor on X faces
    FS5D stressy;       // Stress tensor on Y faces
    FS5D stressz;       // Stress tensor on Z faces
    FS4D gradRho;       // Density Gradient array
    FS4D cFlux;
    FS6D mFlux;
    FS1D cd;            // Device configuration array

};
