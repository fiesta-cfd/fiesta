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

#include "fiesta.hpp"
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
    Fiesta::Simulation sim;
    Fiesta::initialize(sim.cf,argc,argv);

    sim.cf.totalTimer.start();
    sim.cf.initTimer.start();
 
    // Choose Module
    if (sim.cf.ndim == 3){
      if (sim.cf.grid > 0)
        sim.f = std::make_unique<gen3d_func>(sim.cf);
      else
        sim.f = std::make_unique<cart3d_func>(sim.cf);
    }else{
      if (sim.cf.grid > 0)
        sim.f = std::make_unique<gen2d_func>(sim.cf);
      else
        sim.f = std::make_unique<cart2d_func>(sim.cf);
    }

    Fiesta::initializeSimulation(sim);

    Log::message("Executing pre-simulation hook");
    sim.f->preSim();

    sim.cf.initTimer.stop();

    Log::message("Beginning Main Time Loop");
    sim.cf.simTimer.start();
    for (int t = sim.cf.tstart; t < sim.cf.tend+1; ++t) {
      Fiesta::checkIO(sim,t);

      if (sim.cf.exitFlag==1){
        exit_value=1;
        break;
      }

      Fiesta::step(sim,t);
    }
    sim.cf.simTimer.stop();
    Log::message("Simulation complete!");
 
    Log::message("Executing post-simulation hook");
    sim.f->postSim();

    sim.cf.totalTimer.stop();
    Fiesta::reportTimers(sim.cf,sim.f);
  }
  Fiesta::finalize();
  return exit_value;
}
