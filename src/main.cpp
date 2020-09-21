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
 
  // Choose Scheme
  rk_func *f;
  if (cf.ndim == 3) {
    if (cf.grid == 1) {
      f = new gen3d_func(cf);
    } else {
      f = new cart3d_func(cf);
    }
  } else {
    if (cf.grid == 1) {
      f = new gen2d_func(cf);
    } else {
      f = new cart2d_func(cf);
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
