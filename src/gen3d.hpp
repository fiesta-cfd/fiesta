#include "Kokkos_Core.hpp"
#include "fiesta.hpp"
#include "input.hpp"
#include "rkfunction.hpp"

class gen3d_func : public rk_func {

public:
  gen3d_func(struct inputConfig &cf, FS1D &mcd);
  void compute();
  void preStep();
  void postStep();
  void preSim();
  void postSim();

  FS5D metrics;
  FS3D p;       // Pressure
  FS3D T;       // Temperature
  FS4D tvel;    // Transformed Velocity
  FS3D rho;     // Total Density
  FS3D qx;      // Heat Fluxes in X direciton
  FS3D qy;      // Heat Fluxes in Y direction
  FS3D qz;      // Heat Fluxes in Z direction
  FS3D fluxx;   // Weno Fluxes in X direction
  FS3D fluxy;   // Weno Fluxes in Y direction
  FS3D fluxz;   // Weno Fluxes in Z direction
  FS5D stressx; // Stress Tensor on X faces
  FS5D stressy; // Stress Tensor on Y faces
  FS5D stressz; // Stress Tensor on Z faces
  FS4D gradRho; // Density Gradient array
  FS4D cFlux;
  FS6D mFlux;
  FS1D cd; // Device configuration array
#ifndef NOMPI
  FS5D ls, lr, rs, rr, bs, br, ts, tr, hs, hr, fs, fr;
  FS5DH lsH, lrH, rsH, rrH, bsH, brH, tsH, trH, hsH, hrH, fsH, frH;
#endif
};
