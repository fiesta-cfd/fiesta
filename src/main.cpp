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

    // Choose Module
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

    class blockWriter<double> myblock(cf, f, cf.autoRestartName, cf.pathName, false, cf.restart_freq,!cf.autoRestart);

    std::vector<blockWriter<float> > testblocks;
    luaReader L(cf.inputFname,"fiesta");
    L.getIOBlock(cf,f,cf.ndim,testblocks);
    L.close();
 
    // Pre Simulation
    Fiesta::Log::message("Executing pre-simulation hook");
    applyBCs(cf, f);
    f->preSim();

    cf.initTimer.stop();
    cf.simTimer.reset();

    // Main time loop
    Fiesta::Log::message("Beginning Main Time Loop");
    for (int t = cf.tstart; t < cf.tend+1; ++t) {
      Fiesta::collectSignals(cf);

      Fiesta::checkIO(cf,f,t,cf.time,testblocks,myblock);

      if (cf.exitFlag==1)
        break;

      f->preStep();
      rkAdvance(cf,f);
      f->postStep();

      cf.time += cf.dt;
      cf.t = t + 1;
    }
    Fiesta::Log::message("Simulation complete!");
 
    // Post Simulation
    Fiesta::Log::message("Executing post-simulation hook");
    f->postSim();


    // Stop Timers and Report
    cf.simTimer.stop();
    cf.totalTimer.stop();
    Fiesta::Log::message("Reporting Timers:");
    Fiesta::reportTimers(cf,f);
  }
  Fiesta::finalize();
  return 0;
}
