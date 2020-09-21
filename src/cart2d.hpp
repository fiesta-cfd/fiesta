#include "fiesta.hpp"
#include "input.hpp"
#include "rkfunction.hpp"
#ifndef NOMPI
#include "mpi.hpp"
#endif
#include "advect.hpp"
#include "bc.hpp"
#include "buoyancy.hpp"
#include "ceq.hpp"
#include "flux.hpp"
#include "noise.hpp"
#include "particle.hpp"
#include "presgrad.hpp"
#include "secondary.hpp"
#include "velocity.hpp"
#include "viscosity.hpp"

class cart2d_func : public rk_func {

public:
  cart2d_func(struct inputConfig &cf_);

  void compute();
  void preStep();
  void postStep();
  void preSim();
  void postSim();

  FS2D p;       // Pressure
  FS3D vel;     // Velocity
  FS2D T;       // Temperature
  FS2D rho;     // Total Density
  FS2D qx;      // Heat Fluxes in X direction
  FS2D qy;      // Heat Fluxes in X direction
  FS2D fluxx;   // Weno Fluxes in X direction
  FS2D fluxy;   // Weno Fluxes in Y direction
  FS4D stressx; // stress tensor on x faces
  FS4D stressy; // stress tensor on y faces
  FS3D gradRho; // Density Gradient array
  FS2D_I noise; // Noise indicator array
  //    FSP2D particles;   // Particle array
  //    FSP2DH particlesH; // Particle host array
};
