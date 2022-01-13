/*
  Copyright 2019-2021 The University of New Mexico

  This file is part of FIESTA.
  
  FIESTA is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option) any
  later version.
  
  FIESTA is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
  details.
  
  You should have received a copy of the GNU Lesser General Public License
  along with FIESTA.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "kokkosTypes.hpp"
#include "input.hpp"
#include "rkfunction.hpp"
#include "fiesta.hpp"
#include "rk.hpp"
#include "bc.hpp"
#include "debug.hpp"
#include "log.hpp"
#include "block.hpp"
#include "signal.hpp"

#include <vector>
#include "luaReader.hpp"
#include "log2.hpp"

// Compute Objects
#include "cart2d.hpp"
#include "cart3d.hpp"
#include "gen2d.hpp"
#include "gen3d.hpp"

int main(int argc, char *argv[]) {
  int exit_value=0;
  {
    // read input file and initialize configuration
    fsconf cf;
    Fiesta::initialize(cf,argc,argv);

    // create signal handler
    class fiestaSignalHandler *signalHandler = 0;
    signalHandler = signalHandler->getInstance(cf);
    signalHandler->registerSignals();

    cf.totalTimer.start();
    cf.initTimer.start();

    Fiesta::Simulation sim;

    // Choose Module
    /* rk_func *f; */
    if (cf.ndim == 3){
      if (cf.grid > 0){
        sim.f = std::make_unique<gen3d_func>(cf);
      }else{
        sim.f = std::make_unique<cart3d_func>(cf);
      }
    }else{
      if (cf.grid > 0){
        sim.f = std::make_unique<gen2d_func>(cf);
      }else{
        sim.f = std::make_unique<cart2d_func>(cf);
      }
    }
 
    // Initialize Simulation
    Fiesta::initializeSimulation(cf,sim.f);

    /* class blockWriter<double> myblock(cf, sim.f, cf.autoRestartName, cf.pathName, false, cf.restart_freq,!cf.autoRestart); */
    sim.restartview = std::make_unique<blockWriter<double>>(cf, sim.f, cf.autoRestartName, cf.pathName, false, cf.restart_freq,!cf.autoRestart);

    /* std::vector<blockWriter<float> > testblocks; */
    luaReader L(cf.inputFname,"fiesta");
    L.getIOBlock(cf,sim.f,cf.ndim,sim.ioviews);
    L.close();
 
    // Pre Simulation
    Log::message("Executing pre-simulation hook");
    applyBCs(cf, sim.f);
    sim.f->preSim();

    cf.initTimer.stop();
    cf.simTimer.reset();

    // Main time loop
    Log::message("Beginning Main Time Loop");
    for (int t = cf.tstart; t < cf.tend+1; ++t) {
      Fiesta::collectSignals(cf);

      Fiesta::checkIO(cf,sim.f,t,cf.time,sim.ioviews,sim.restartview);

      if (cf.exitFlag==1){
        exit_value=1;
        break;
      }

      sim.f->preStep();
      rkAdvance(cf,sim.f);
      sim.f->postStep();

      cf.time += cf.dt;
      cf.t = t + 1;
    }
    Log::message("Simulation complete!");
 
    // Post Simulation
    Log::message("Executing post-simulation hook");
    sim.f->postSim();


    // Stop Timers and Report
    cf.simTimer.stop();
    cf.totalTimer.stop();
    Log::message("Reporting Timers:");
    Fiesta::reportTimers(cf,sim.f);
  }
  Fiesta::finalize();
  return exit_value;
}
