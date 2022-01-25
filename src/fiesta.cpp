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
#ifdef HAVE_MPI
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
#include "reader.hpp"
#include "yogrt.h"
#include "signal.hpp"
#include "rk.hpp"

using namespace std;

//
// Initialize Fiesta and fill configuration structure with input file variables
//
struct inputConfig Fiesta::initialize(struct inputConfig &cf, int argc, char **argv){
  struct commandArgs cArgs = getCommandlineOptions(argc, argv);

  // Initialize MPI and get temporary rank.
  int temp_rank = 0;
#ifdef HAVE_MPI
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &temp_rank);
#endif

  if (temp_rank == 0) printSplash(cArgs.colorFlag);

  Log::Logger(cArgs.verbosity,cArgs.colorFlag,temp_rank);

  // Initialize Kokkos
  Log::message("Initializing Kokkos");
  Kokkos::InitArguments kokkosArgs;
#ifdef HAVE_CUDA
  kokkosArgs.ndevices = cArgs.numDevices;
#elif HAVE_OPENMP
  kokkosArgs.num_threads = cArgs.numThreads;
#endif
  Kokkos::initialize(kokkosArgs);

  // Execute lua script and get input parameters
  Log::message("Executing Lua Input Script");
  executeConfiguration(cf,cArgs);

#ifdef HAVE_MPI
  // perform domain decomposition
  Log::message("Initializing MPI Setup");
  mpi_init(cf);
  MPI_Barrier(cf.comm);
#endif
  Log::message("Printing Configuration");
  printConfig(cf);

  // create signal handler
  class fiestaSignalHandler *signalHandler = 0;
  signalHandler = signalHandler->getInstance(cf);
  signalHandler->registerSignals();

  return cf;
}

//
// Initialize the simulation and load initial data
//
/* void Fiesta::initializeSimulation(struct inputConfig &cf, std::unique_ptr<class rk_func>&f){ */
void Fiesta::initializeSimulation(Simulation &sim){
#ifdef HAVE_MPI
  if (sim.cf.visc || sim.cf.ceq){
    if (sim.cf.mpiScheme == 1)
      sim.cf.m = std::make_shared<orderedHaloExchange>(sim.cf,sim.f->var);
    else if (sim.cf.mpiScheme == 2)
      sim.cf.m = std::make_shared<orderedHostHaloExchange>(sim.cf,sim.f->var);
    else{
      Log::error("Invalid MPI type.  Only 'gpu-aware' and 'host' are available when viscosity or c-equations are enabled.");
      exit(EXIT_FAILURE);
    }
  }else{
    if (sim.cf.mpiScheme == 1)
      sim.cf.m = std::make_shared<copyHaloExchange>(sim.cf,sim.f->var);
    else if (sim.cf.mpiScheme == 2)
      sim.cf.m = std::make_shared<packedHaloExchange>(sim.cf,sim.f->var);
    else if (sim.cf.mpiScheme == 3)
      sim.cf.m = std::make_shared<directHaloExchange>(sim.cf,sim.f->var);
  }
#endif
  //sim.cf.w = std::make_shared<hdfWriter>(sim.cf,f);


  // If not restarting, generate initial conditions and grid
  if (sim.cf.restart == 0) {
    // Generate Grid Coordinates
    Log::message("Generating grid");
    sim.cf.gridTimer.start();
    loadGrid(sim.cf, sim.f->grid);
    sim.cf.gridTimer.stop();
    Log::message("Grid generated in: {}",sim.cf.gridTimer.get());

    // Generate Initial Conditions
    Log::message("Generating initial conditions");
    sim.cf.loadTimer.start();
    loadInitialConditions(sim.cf, sim.f->var, sim.f->grid);
    sim.cf.loadTimer.stop();
    Log::message("Initial conditions generated in: {}",sim.cf.loadTimer.get());

    // sim.cf.writeTimer.start();
    // // Write initial solution file
    // if (sim.cf.write_freq > 0) {
    //   sim.f->timers["solWrite"].reset();
    //   sim.cf.w->writeSolution(sim.cf, f, 0, 0.00);
    //   sim.f->timers["solWrite"].accumulate();
    // }

    // // Write Initial Restart File
    // if (sim.cf.restart_freq > 0) {
    //   sim.f->timers["resWrite"].reset();
    //   sim.cf.w->writeRestart(sim.cf, f, 0, 0.00);
    //   sim.f->timers["resWrite"].accumulate();
    // }
    // sim.cf.writeTimer.stop();
  }else{ // If Restarting, Load Restart File
    sim.cf.writeTimer.start();
    Log::message("Loading restart file:");
    sim.cf.loadTimer.reset();
    readRestart(sim.cf, sim.f);
    sim.cf.loadTimer.stop();
    Log::message("Loaded restart data in: {}",sim.cf.loadTimer.get());
  }

  if (sim.cf.rank==0){
    if (!std::filesystem::exists(sim.cf.pathName)){
      Log::message("Creating directory: '{}'",sim.cf.pathName);
      std::filesystem::create_directories(sim.cf.pathName);
    }
  }

  sim.restartview = std::make_unique<blockWriter<FSCAL>>(sim.cf, sim.f, sim.cf.autoRestartName, sim.cf.pathName, false, sim.cf.restart_freq,!sim.cf.autoRestart);

  luaReader L(sim.cf.inputFname,"fiesta");
  L.getIOBlock(sim.cf,sim.f,sim.cf.ndim,sim.ioviews);
  L.close();

}

// Write solutions, restarts and status checks
void Fiesta::checkIO(Simulation &sim, size_t t){
  collectSignals(sim.cf);
  // Print current time step
  if (sim.cf.rank == 0) {
    if (sim.cf.out_freq > 0)
      if (t % sim.cf.out_freq == 0)
        Log::info("[{}] Timestep {} of {}. Simulation Time: {:.2e}s",t,t,sim.cf.tend,sim.cf.time);
  }

  // Print status check if necessary
  if (sim.cf.stat_freq > 0) {
    if (t % sim.cf.stat_freq == 0) {
      sim.f->timers["statCheck"].reset();
      statusCheck(sim.cf.colorFlag, sim.cf, sim.f, sim.cf.time, sim.cf.totalTimer, sim.cf.simTimer);
      sim.f->timers["statCheck"].accumulate();
    }
  }

  // Write solution file if necessary
  if (sim.cf.write_freq > 0) {
    if (t % sim.cf.write_freq == 0) {
      sim.f->timers["solWrite"].reset();
      sim.cf.w->writeSolution(sim.cf, sim.f, t, sim.cf.time);
      Kokkos::fence();
      sim.f->timers["solWrite"].accumulate();
    }
  }

  // Check Restart Frequency
  if (sim.cf.restart_freq > 0 && t > sim.cf.tstart) {
    if (t % sim.cf.restart_freq == 0) {
      sim.cf.restartFlag=1;
    }
  }

  // Check Time-Remaining and Flag Restart
  if (sim.cf.rank==0){
    if(yogrt_remaining() < sim.cf.restartTimeRemaining){
      sim.cf.restartFlag = 1;
      sim.cf.exitFlag = 1;
      Log::error("Time remaining is less than {}s:  Writing restart and exiting after timestep {}.",sim.cf.restartTimeRemaining,sim.cf.t);
    }
  }

  #ifdef HAVE_MPI
  {
  int glblRestartFlag=0;
  int glblExitFlag=0;
  MPI_Allreduce(&sim.cf.restartFlag,&glblRestartFlag,1,MPI_INT,MPI_MAX,sim.cf.comm);
  sim.cf.restartFlag=glblRestartFlag;
  MPI_Allreduce(&sim.cf.exitFlag,&glblExitFlag,1,MPI_INT,MPI_MAX,sim.cf.comm);
  sim.cf.exitFlag=glblExitFlag;
  }
  #endif

  // Write restart file if necessary
  if (sim.cf.restartFlag==1){
    sim.f->timers["resWrite"].reset();
    //sim.cf.w->writeRestart(sim.cf, f, t, time);
    sim.restartview->write(sim.cf,sim.f,t,sim.cf.time);
    Kokkos::fence();
    sim.f->timers["resWrite"].accumulate();
    sim.cf.restartFlag=0;
  }

  // Write solution blocks
  for (auto& block : sim.ioviews){
    if(block.frq() > 0){
      if (t % block.frq() == 0) {
        sim.f->timers["solWrite"].reset();
        block.write(sim.cf,sim.f,t,sim.cf.time);
        sim.f->timers["solWrite"].accumulate();
      }
    }
    sim.cf.ioThisStep = false;
    if (t % block.frq() == 0){
      sim.cf.ioThisStep = true;
    }
  }
}

void Fiesta::collectSignals(struct inputConfig &cf){
#ifdef HAVE_MPI
  int glblRestartFlag=0;
  int glblExitFlag=0;

  MPI_Allreduce(&cf.restartFlag,&glblRestartFlag,1,MPI_INT,MPI_MAX,cf.comm);
  cf.restartFlag=glblRestartFlag;

  MPI_Allreduce(&cf.exitFlag,&glblExitFlag,1,MPI_INT,MPI_MAX,cf.comm);
  cf.exitFlag=glblExitFlag;
#endif

  if (cf.restartFlag && cf.exitFlag)
      Log::warning("Recieved SIGURG:  Writing restart and exiting after timestep {}.",cf.t);
  else if (cf.restartFlag)
      Log::message("Recieved SIGUSR1:  Writing restart after timestep {}.",cf.t);
  else if (cf.exitFlag)
      Log::error("Recieved SIGTERM:  Exiting after timestep {}.",cf.t);
}

void Fiesta::reportTimers(struct inputConfig &cf, std::unique_ptr<class rk_func>&f){
  Log::message("Reporting Timers:");
  // Sort computer timers
  typedef std::function<bool(std::pair<std::string, Timer::fiestaTimer>,
                             std::pair<std::string, Timer::fiestaTimer>)>
      Comparator;
  Comparator compFunctor = [](std::pair<std::string, Timer::fiestaTimer> elem1,
                              std::pair<std::string, Timer::fiestaTimer> elem2) {
    return elem1.second.get() > elem2.second.get();
  };
  std::set<std::pair<std::string, Timer::fiestaTimer>, Comparator> stmr(
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

void Fiesta::step(Simulation &sim, size_t t){
      sim.f->preStep();
      rkAdvance(sim.cf,sim.f);
      sim.f->postStep();

      sim.cf.time += sim.cf.dt;
      sim.cf.t = t + 1;
}
 
// clean up kokkos and mpi
void Fiesta::finalize(){
  Kokkos::finalize();
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
