#include "fiesta.hpp"
#include "input.hpp"
#include "noise.hpp"
#include "rkfunction.hpp"
#ifndef NOMPI
#include "mpi.hpp"
#endif
#include "advect.hpp"
#include "flux.hpp"
#include "metric.hpp"
#include "particle.hpp"
#include "presgrad.hpp"
#include "secondary.hpp"
#include "velocity.hpp"

class gen2d_func : public rk_func {

public:
  gen2d_func(struct inputConfig &cf_, FS1D &cd_);

  void compute();
  void preStep();
  void postStep();
  void preSim();
  void postSim();

  FS2D p;       // Pressure
  FS2D T;       // Temperature
  FS3D tvel;    // Transformed velocity
  FS2D rho;     // Total Density
  FS2D fluxx;   // Weno Fluxes in X direction
  FS2D fluxy;   // Weno Fluxes in Y direction
  FS2D_I noise; // Noise indicator array
  FS1D cd;      // Device configuration array
  FS4D metrics; // jacobian metrics
  //    FSP2D particles;   // Particle array
  //    FSP2DH particlesH; // Host particle array

#ifndef NOMPI
  FS4D ls, lr, rs, rr, bs, br, ts, tr;
  FS4DH lsH, lrH, rsH, rrH, bsH, brH, tsH, trH;
#endif
};
