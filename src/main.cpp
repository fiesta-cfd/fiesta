#include "kokkosTypes.hpp"
#include "input.hpp"
#include "rkfunction.hpp"
#include "fiesta.hpp"
#include "rk.hpp"
#include "debug.hpp"

// Compute Objects
#include "cart2d.hpp"
#include "cart3d.hpp"
#include "gen2d.hpp"
#include "gen3d.hpp"

int main(int argc, char *argv[]) {
  
  struct inputConfig cf = Fiesta::initialize(argc,argv);

  cf.totalTimer.start();
  cf.initTimer.start();
 
  // Choose Scheme
  rk_func *f;
  if (cf.ndim == 3){
    if (cf.grid > 0){
      f = new gen3d_func(cf);
    }else{
      f = new cart3d_func(cf);
    }
  }else{
    if (cf.grid > 0){
      f = new gen2d_func(cf);
    }else{
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

  // Stop Timers and Report
  cf.simTimer.stop();
  cf.totalTimer.stop();
  Fiesta::reportTimers(cf,f);
 
  Fiesta::finalize();
  return 0;
}
