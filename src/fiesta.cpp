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

#include <memory>
#include "fiesta.hpp"
#include "input.hpp"
#include "debug.hpp"
#ifndef NOMPI
#include "mpi.hpp"
#endif
#include "output.hpp"
#include "rkfunction.hpp"
#include "status.hpp"
#include <set>
#include "log.hpp"
#include "block.hpp"
//#include <csignal>
#include <string>
#include <vector>
#include "luaReader.hpp"
#include "fmt/core.h"
#include "pretty.hpp"
#include <filesystem>
#include "log2.hpp"
#include <iostream>

using namespace std;

int Fiesta::fiestaTest(int a, int b){
    return a + b;
}

//
// Initialize Fiesta and fill configuration structure with input file variables
//
struct inputConfig Fiesta::initialize(struct inputConfig &cf, int argc, char **argv){
  struct commandArgs cArgs = getCommandlineOptions(argc, argv);

  // Initialize MPI and get temporary rank.
  int temp_rank = 0;
#ifndef NOMPI
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &temp_rank);
#endif

  if (temp_rank == 0) printSplash(cArgs.colorFlag);

  Fiesta::Log::Logger(cArgs.verbosity,cArgs.colorFlag,temp_rank);

  // Initialize Kokkos
  Fiesta::Log::message("Initializing Kokkos");
  Kokkos::InitArguments kokkosArgs;
#ifdef HAVE_CUDA
  kokkosArgs.ndevices = cArgs.numDevices;
#elif HAVE_OPENMP
  kokkosArgs.num_threads = cArgs.numThreads;
#endif
  Kokkos::initialize(kokkosArgs);

  // Execute lua script and get input parameters
  Fiesta::Log::message("Executing Lua Input Script");
  executeConfiguration(cf,cArgs);

#ifndef NOMPI
  // perform domain decomposition
  Fiesta::Log::message("Initializing MPI Setup");
  mpi_init(cf);
  Fiesta::Log::debugAll("Subdomain Dimensions=({},{},{}), Offset=({},{},{})",cf.nci,cf.ncj,cf.nck,cf.subdomainOffset[0],cf.subdomainOffset[1],cf.subdomainOffset[2]);
  MPI_Barrier(cf.comm);
#endif
  Fiesta::Log::message("Printing Configuration");
  printConfig(cf);

  return cf;
}

//
// Initialize the simulation and load initial data
//
void Fiesta::initializeSimulation(struct inputConfig &cf, rk_func *f){
  //create IO object
#ifdef NOMPI
  cf.w = new serialVTKWriter(cf, f->grid, f->var);
#else
  if (cf.mpiScheme == 1)
    cf.m = new copyHaloExchange(cf, f->var);
  else if (cf.mpiScheme == 2)
    cf.m = new packedHaloExchange(cf, f->var);
  else if (cf.mpiScheme == 3)
    cf.m = new directHaloExchange(cf, f->var);
  //cf.w = new hdfWriter(cf, f);
  //cf.m = new mpiBuffers(cf);
  cf.w = std::make_shared<hdfWriter>(cf,f);
  cf.m = std::make_shared<mpiBuffers>(cf);

  //luaReader L(cf.inputFname);
  //L.getIOBlock(cf,f,cf.ndim,cf.ioblocks);
  //L.close();

#endif

  // If not restarting, generate initial conditions and grid
  if (cf.restart == 0) {
    // Generate Grid Coordinates
    Fiesta::Log::message("Generating grid");
    cf.gridTimer.start();
    loadGrid(cf, f->grid);
    cf.gridTimer.stop();
    Fiesta::Log::message("Grid generated in: {}",cf.gridTimer.get());

    // Generate Initial Conditions
    Fiesta::Log::message("Generating initial conditions");
    cf.loadTimer.start();
    loadInitialConditions(cf, f->var, f->grid);
    cf.loadTimer.stop();
    Fiesta::Log::message("Initial conditions generated in: {}",cf.loadTimer.get());

    // cf.writeTimer.start();
    // // Write initial solution file
    // if (cf.write_freq > 0) {
    //   f->timers["solWrite"].reset();
    //   cf.w->writeSolution(cf, f, 0, 0.00);
    //   f->timers["solWrite"].accumulate();
    // }

    // // Write Initial Restart File
    // if (cf.restart_freq > 0) {
    //   f->timers["resWrite"].reset();
    //   cf.w->writeRestart(cf, f, 0, 0.00);
    //   f->timers["resWrite"].accumulate();
    // }
    // cf.writeTimer.stop();
  }else{ // If Restarting, Load Restart File
    cf.writeTimer.start();
    Fiesta::Log::message("Loading restart file:");
    cf.loadTimer.reset();
    cf.w->readSolution(cf, f->grid, f->var);
    cf.loadTimer.stop();
    Fiesta::Log::message("Loaded restart data in: {}",cf.loadTimer.get());
  }

  if (cf.rank==0){
    if (!std::filesystem::exists(cf.pathName)){
      Fiesta::Log::message("Creating directory: '{}'",cf.pathName);
      std::filesystem::create_directories(cf.pathName);
    }
  }
}

// Write solutions, restarts and status checks
void Fiesta::checkIO(struct inputConfig &cf, rk_func *f, int t, double time,vector<blockWriter<float> >& ioblocks, blockWriter<double>& rsblock){
  // Print current time step
  if (cf.rank == 0) {
    if (cf.out_freq > 0)
      if (t % cf.out_freq == 0)
        Fiesta::Log::info("[{}] Timestep {} of {}. Simulation Time: {:.2e}s",t,t,cf.tend,time);
  }

  // Print status check if necessary
  if (cf.stat_freq > 0) {
    if (t % cf.stat_freq == 0) {
      f->timers["statCheck"].reset();
      statusCheck(cf.colorFlag, cf, f, time, cf.totalTimer, cf.simTimer);
      f->timers["statCheck"].accumulate();
    }
  }

  // Write solution file if necessary
  if (cf.write_freq > 0) {
    if (t % cf.write_freq == 0) {
      f->timers["solWrite"].reset();
      cf.w->writeSolution(cf, f, t, time);
      Kokkos::fence();
      f->timers["solWrite"].accumulate();
    }
  }

  // Check Restart Frequency
  if (cf.restart_freq > 0 && t > 0) {
    if (t % cf.restart_freq == 0) {
      cf.restartFlag=1;
    }
  }

  // Write restart file if necessary
  if (cf.restartFlag==1){
    f->timers["resWrite"].reset();
    //cf.w->writeRestart(cf, f, t, time);
    rsblock.write(cf,f,t,time);
    Kokkos::fence();
    f->timers["resWrite"].accumulate();
    cf.restartFlag=0;
  }

  // Write solution blocks
  for (auto& block : ioblocks){
    if(block.frq() > 0){
      if (t % block.frq() == 0) {
        f->timers["solWrite"].reset();
        block.write(cf,f,t,time);
        f->timers["solWrite"].accumulate();
      }
    }
  }
}

void Fiesta::collectSignals(struct inputConfig &cf){
  int glblRestartFlag=0;
  int glblExitFlag=0;

  MPI_Allreduce(&cf.restartFlag,&glblRestartFlag,1,MPI_INT,MPI_MAX,cf.comm);
  cf.restartFlag=glblRestartFlag;

  MPI_Allreduce(&cf.exitFlag,&glblExitFlag,1,MPI_INT,MPI_MAX,cf.comm);
  cf.exitFlag=glblExitFlag;

  if (cf.restartFlag && cf.exitFlag)
      Fiesta::Log::warning("Recieved SIGURG:  Writing restart and exiting after timestep {}.",cf.t);
  else if (cf.restartFlag)
      Fiesta::Log::message("Recieved SIGUSR1:  Writing restart after timestep {}.",cf.t);
  else if (cf.exitFlag)
      Fiesta::Log::error("Recieved SIGTERM:  Exiting after timestep {}.",cf.t);
}

void Fiesta::reportTimers(struct inputConfig &cf, rk_func *f){
  // Sort computer timers
  typedef std::function<bool(std::pair<std::string, fiestaTimer>,
                             std::pair<std::string, fiestaTimer>)>
      Comparator;
  Comparator compFunctor = [](std::pair<std::string, fiestaTimer> elem1,
                              std::pair<std::string, fiestaTimer> elem2) {
    return elem1.second.get() > elem2.second.get();
  };
  std::set<std::pair<std::string, fiestaTimer>, Comparator> stmr(
      f->timers.begin(), f->timers.end(), compFunctor);

  if (cf.rank==0){
    using fmt::format;
    ansiColors c(cf.colorFlag);

    string timerFormat = format("{: >8}{{}}{{: <{}}}{{:{}.3e}}{}\n","",32,16,c(reset));
    cout << format(timerFormat,c(green),"Total Execution Time",cf.totalTimer.get());
    cout << "\n";

    cout << format(timerFormat,c(green),"Total Startup Time",cf.initTimer.get());
    if (cf.restart==1)
      cout << format(timerFormat,c(reset),"Restart Read",cf.loadTimer.get());
    else{
      cout << format(timerFormat,c(reset),"Initial Condition Generation",cf.loadTimer.get());
      cout << format(timerFormat,c(reset),"Grid Generation",cf.gridTimer.get());
      cout << format(timerFormat,c(reset),"Initial Condition Write Time",cf.writeTimer.get());
    }
    cout << "\n";

    cout << format(timerFormat,c(green),"Total Simulation Time",cf.simTimer.get());
    for (auto tmr : stmr)
      cout << format(timerFormat,c(reset),tmr.second.describe(),tmr.second.get());
  }
}
 
// clean up kokkos and mpi
void Fiesta::finalize(){
  Kokkos::finalize();
#ifndef NOMPI
  MPI_Finalize();
#endif
}
