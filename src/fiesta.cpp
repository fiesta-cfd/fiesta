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
#include "input.hpp"
#include "debug.hpp"
#ifndef NOMPI
#include "mpi.hpp"
#endif
#include "output.hpp"
#include "rkfunction.hpp"
#include "status.hpp"
#include <iostream>
#include <set>

using namespace std;

//
// Initialize Fiesta and fill configuration structure with input file variables
//
struct inputConfig Fiesta::initialize(int argc, char **argv){
  struct commandArgs cArgs = getCommandlineOptions(argc, argv);

  Kokkos::InitArguments kokkosArgs;
#ifdef HAVE_CUDA
  kokkosArgs.ndevices = cArgs.numDevices;
#elif HAVE_OPENMP
  kokkosArgs.num_threads = cArgs.numThreads;
#endif

  // Initialize MPI and get temporary rank.
  int temp_rank = 0;
#ifndef NOMPI
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &temp_rank);
#endif

  if (temp_rank == 0) printSplash(cArgs.colorFlag);

  // Initialize kokkos and set kokkos finalize as exit function.
  Kokkos::initialize(kokkosArgs);
  // Kokkos::initialize(argc, argv);

  // Execute lua script and get input parameters
  struct inputConfig cf = executeConfiguration(cArgs);

#ifndef NOMPI
  // perform domain decomposition and MPI initialization.
  mpi_init(cf);
#endif
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
  cf.w = new hdfWriter(cf, f);
  if (cf.mpiScheme == 1)
    cf.m = new copyHaloExchange(cf, f->var);
  else if (cf.mpiScheme == 2)
    cf.m = new packedHaloExchange(cf, f->var);
  else if (cf.mpiScheme == 3)
    cf.m = new directHaloExchange(cf, f->var);
#endif

  // If not restarting, generate initial conditions and grid
  if (cf.restart == 0) {
    // Generate Grid Coordinates
    if (cf.rank == 0)
      cout << c(cf.colorFlag, GRE) << "Generating Grid:" << c(cf.colorFlag, NON) << endl;
    cf.gridTimer.start();
    loadGrid(cf, f->grid);
    cf.gridTimer.stop();
    if (cf.rank == 0)
      cout << "    Generated in: " << c(cf.colorFlag, CYA) << cf.gridTimer.getf(cf.timeFormat)
           << c(cf.colorFlag, NON) << endl
           << endl;

    // Generate Initial Conditions
    if (cf.rank == 0)
      cout << c(cf.colorFlag, GRE)
           << "Generating Initial Conditions:" << c(cf.colorFlag, NON) << endl;
    cf.loadTimer.start();
    loadInitialConditions(cf, f->var, f->grid);
    cf.loadTimer.stop();
    if (cf.rank == 0)
      cout << "    Generated in: " << c(cf.colorFlag, CYA) << cf.loadTimer.getf(cf.timeFormat)
           << c(cf.colorFlag, NON) << endl
           << endl;

    cf.writeTimer.start();
    // Write Initial Solution File
    if (cf.rank == 0)
      if (cf.write_freq > 0 || cf.restart_freq > 0)
        cout << c(cf.colorFlag, GRE)
             << "Writing Initial Conditions:" << c(cf.colorFlag, NON) << endl;
    if (cf.write_freq > 0) {
      f->timers["solWrite"].reset();
      cf.w->writeSolution(cf, f, 0, 0.00);
      f->timers["solWrite"].accumulate();
    }
    // Write Initial Restart File
    if (cf.restart_freq > 0) {
      f->timers["resWrite"].reset();
      cf.w->writeRestart(cf, f, 0, 0.00);
      f->timers["resWrite"].accumulate();
    }
    cf.writeTimer.stop();
    if (cf.rank == 0)
      if (cf.write_freq > 0 || cf.restart_freq > 0)
        cout << "    Wrote in: " << c(cf.colorFlag, CYA) << cf.writeTimer.getf(cf.timeFormat)
             << c(cf.colorFlag, NON) << endl;

  }else{ // If Restarting, Load Restart File
    cf.writeTimer.start();
    if (cf.rank == 0)
      cout << c(cf.colorFlag, GRE) << "Loading Restart File:" << c(cf.colorFlag, NON)
           << endl;
    cf.loadTimer.reset();
    cf.w->readSolution(cf, f->grid, f->var);
    cf.loadTimer.stop();
    if (cf.rank == 0)
      cout << "    Loaded in: " << setprecision(2) << c(cf.colorFlag, CYA)
           << cf.loadTimer.getf(cf.timeFormat) << "s" << c(cf.colorFlag, NON) << endl;

  }
  // notify simulation start
  if (cf.rank == 0) {
    cout << endl << "-----------------------" << endl << endl;
   cout << c(cf.colorFlag, GRE) << "Starting Simulation:" << c(cf.colorFlag, NON) << endl;
  }
}

// Write solutions, restarts and status checks
void Fiesta::checkIO(struct inputConfig &cf, rk_func *f, int t, double time){
  // Print current time step
  if (cf.rank == 0) {
    if (cf.out_freq > 0)
      if ((t + 1) % cf.out_freq == 0)
        cout << c(cf.colorFlag, YEL) << left << setw(15)
             << "    Iteration:" << c(cf.colorFlag, NON) << c(cf.colorFlag, CYA) << right
             << setw(0) << t + 1 << c(cf.colorFlag, NON) << "/" << c(cf.colorFlag, CYA)
             << left << setw(0) << cf.tend << c(cf.colorFlag, NON) << ", "
             << c(cf.colorFlag, CYA) << right << setw(0) << setprecision(3)
             << scientific << time << "s" << c(cf.colorFlag, NON) << endl;
  }
  // Write solution file if necessary
  if (cf.write_freq > 0) {
    if ((t + 1) % cf.write_freq == 0) {
      f->timers["solWrite"].reset();
      cf.w->writeSolution(cf, f, t + 1, time);
      Kokkos::fence();
      f->timers["solWrite"].accumulate();
    }
  }
  // Write restart file if necessary
  if (cf.restart_freq > 0) {
    if ((t + 1) % cf.restart_freq == 0) {
      f->timers["resWrite"].reset();
      cf.w->writeRestart(cf, f, t + 1, time);
      Kokkos::fence();
      f->timers["resWrite"].accumulate();
    }
  }
  // Print status check if necessary
  if (cf.stat_freq > 0) {
    if ((t + 1) % cf.stat_freq == 0) {
      f->timers["statCheck"].reset();
      statusCheck(cf.colorFlag, cf, f, time, cf.totalTimer, cf.simTimer);
      f->timers["statCheck"].accumulate();
    }
  }
}

void Fiesta::reportTimers(struct inputConfig &cf, rk_func *f){
  // notify simulation complete
  if (cf.rank == 0)
    cout << c(cf.colorFlag, GRE) << "Simulation Complete!" << c(cf.colorFlag, NON) << endl;

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

  // print timer values
  if (cf.rank == 0) {
    cout << endl << "-----------------------" << endl << endl;
    cout.precision(2);
    cout << c(cf.colorFlag, GRE) << left << setw(36)
         << "Total Time:" << c(cf.colorFlag, NON) << c(cf.colorFlag, CYA) << right
         << setw(13) << cf.totalTimer.getf(cf.timeFormat) << c(cf.colorFlag, NON) << endl
         << endl;

    cout << c(cf.colorFlag, GRE) << left << setw(36)
         << "  Setup Time:" << c(cf.colorFlag, NON) << c(cf.colorFlag, CYA) << right
         << setw(13) << cf.initTimer.getf(cf.timeFormat) << c(cf.colorFlag, NON) << endl;
    if (cf.restart == 1)
      cout << c(cf.colorFlag, NON) << left << setw(36)
           << "    Restart Read:" << c(cf.colorFlag, NON) << c(cf.colorFlag, CYA) << right
           << setw(13) << cf.loadTimer.getf(cf.timeFormat) << c(cf.colorFlag, NON) << endl
           << endl;
    else {
      cout << c(cf.colorFlag, NON) << left << setw(36)
           << "    Initial Condition Generation:" << c(cf.colorFlag, NON)
           << c(cf.colorFlag, CYA) << right << setw(13) << cf.loadTimer.getf(cf.timeFormat)
           << c(cf.colorFlag, NON) << endl;
      cout << c(cf.colorFlag, NON) << left << setw(36)
           << "    Grid Generation:" << c(cf.colorFlag, NON) << c(cf.colorFlag, CYA)
           << right << setw(13) << cf.gridTimer.getf(cf.timeFormat) << c(cf.colorFlag, NON) << endl;
      cout << c(cf.colorFlag, NON) << left << setw(36)
           << "    Initial Condition WriteTime:" << c(cf.colorFlag, NON)
           << c(cf.colorFlag, CYA) << right << setw(13) << cf.writeTimer.getf(cf.timeFormat)
           << c(cf.colorFlag, NON) << endl
           << endl;
    }

    cout << c(cf.colorFlag, GRE) << left << setw(36)
         << "  Simulation Time:" << c(cf.colorFlag, NON) << c(cf.colorFlag, CYA) << right
         << setw(13) << cf.simTimer.getf(cf.timeFormat) << c(cf.colorFlag, NON) << endl;
    for (auto tmr : stmr) {
      cout << c(cf.colorFlag, NON) << left << setw(36)
           << "    " + tmr.second.describe() + ":" << c(cf.colorFlag, NON)
           << c(cf.colorFlag, CYA) << right << setw(13) << tmr.second.getf(cf.timeFormat)
           << c(cf.colorFlag, NON) << endl;
    }
    cout << " " << endl;
  }
}
 
// clean up kokkos and mpi
void Fiesta::finalize(){
  Kokkos::finalize();
#ifndef NOMPI
  MPI_Finalize();
#endif
}
