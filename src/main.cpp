#include "fiesta.hpp"
#include "input.hpp"
#ifndef NOMPI
#include "hdf.hpp"
#include "mpi.hpp"
#include "mpi.h"
#else
#include "vtk.hpp"
#endif
#include "Kokkos_Core.hpp"
#include "bc.hpp"
#include "cart2d.hpp"
#include "cart3d.hpp"
#include "debug.hpp"
#include "gen2d.hpp"
#include "gen3d.hpp"
#include "output.hpp"
#include "rkfunction.hpp"
#include "status.hpp"
#include "timer.hpp"
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <set>
#include "luaReader.hpp"
#include "fiesta2.hpp"
#include "rk.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  

  struct inputConfig cf = Fiesta::initialize(argc,argv);

  cf.totalTimer.start();
  cf.initTimer.start();

  // Create and copy minimal configuration array for data needed
  // withing Kokkos kernels.
  FS1D cd("deviceCF", 6 + cf.ns * 3);
  typename Kokkos::View<double *>::HostMirror hostcd =
      Kokkos::create_mirror_view(cd);
  Kokkos::deep_copy(hostcd, cd);
  hostcd(0) = cf.ns; // number of gas species
  hostcd(1) = cf.dx; // cell size
  hostcd(2) = cf.dy; // cell size
  hostcd(3) = cf.dz; // cell size
  hostcd(4) = cf.nv; // number of flow variables
  hostcd(5) = cf.ng; // number of flow variables
 
  // include gas properties for each gas species
  int sdx = 6;
  for (int s = 0; s < cf.ns; ++s) {
    hostcd(sdx) = cf.gamma[s];        // ratio of specific heats
    hostcd(sdx + 1) = cf.R / cf.M[s]; // species gas comstant
    hostcd(sdx + 2) = cf.mu[s];       // kinematic viscosity
    sdx += 3;
  }
  Kokkos::deep_copy(cd, hostcd); // copy congifuration array to device
 
  // Choose Scheme
  rk_func *f;
  if (cf.ndim == 3) {
    if (cf.grid == 1) {
      f = new gen3d_func(cf, cd);
    } else {
      f = new cart3d_func(cf,cd);
    }
  } else {
    if (cf.grid == 1) {
      f = new gen2d_func(cf, cd);
    } else {
      f = new cart2d_func(cf, cd);
    }
  }
 
  // Initialize Simulation
  Fiesta::initializeSimulation(cf,f);
 
  // Pre Simulation
  f->preSim();
  cf.initTimer.stop();
  cf.simTimer.reset();
 
  // Main time loop
  for (int t = cf.tstart; t < cf.tend; ++t) {
    cf.time += cf.dt;
    cf.t = t + 1;
 
    f->preStep();
    rkAdvance(cf,f);
    f->postStep();
 
    Fiesta::checkIO(cf,f,t,cf.time);
  }
 
  // Post Simulation
  f->postSim();
  cf.simTimer.stop();
  Fiesta::reportTimers(cf,f);
 
  Fiesta::finalize();
  return 0;
}
