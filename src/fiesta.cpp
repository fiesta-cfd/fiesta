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

struct inputConfig Fiesta::initialize(int argc, char **argv){
  struct commandArgs cArgs = getCommandlineOptions(argc, argv);

  // Initialize MPI and get temporary rank.
  int temp_rank = 0;
#ifndef NOMPI
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &temp_rank);
#endif

  if (temp_rank == 0) printSplash(cArgs.colorFlag);

  // Initialize kokkos and set kokkos finalize as exit function.
  Kokkos::initialize(argc, argv);

  // Execute input file and generate simulation configuration
  struct inputConfig cf;
  cf = executeConfiguration(cArgs);

#ifndef NOMPI
  mpi_init(cf);
#endif
  printConfig(cf);

  return cf;
}

void Fiesta::initializeSimulation(struct inputConfig &cf, rk_func *f){
#ifdef NOMPI
    cf.w = new serialVTKWriter(cf, f->grid, f->var);
#else
    cf.w = new hdfWriter(cf, f);
    cf.m = new mpiBuffers(cf);
#endif

  // If not restarting, generate initial conditions and grid
  if (cf.restart == 0) {
    if (cf.rank == 0)
      cout << c(cf.colorFlag, GRE)
           << "Generating Initial Conditions:" << c(cf.colorFlag, NON) << endl;
    cf.loadTimer.start();
    loadInitialConditions(cf, f->var);
    cf.loadTimer.stop();
    if (cf.rank == 0)
      cout << "    Generated in: " << c(cf.colorFlag, CYA) << cf.loadTimer.getf(cf.timeFormat)
           << c(cf.colorFlag, NON) << endl
           << endl;

    if (cf.rank == 0)
      cout << c(cf.colorFlag, GRE) << "Generating Grid:" << c(cf.colorFlag, NON) << endl;
    f->timers["gridTimer"].start();
    loadGrid(cf, f->grid);
    f->timers["gridTimer"].stop();
    if (cf.rank == 0)
      cout << "    Generated in: " << c(cf.colorFlag, CYA) << f->timers["gridTimer"].getf(cf.timeFormat)
           << c(cf.colorFlag, NON) << endl
           << endl;

    if (cf.particle == 1) {
      if (cf.rank == 0)
        cout << c(cf.colorFlag, GRE) << "Generating Particles:" << c(cf.colorFlag, NON)
             << endl;
      f->timers["gridTimer"].start();
      loadParticles(cf, f->particles);
      f->timers["gridTimer"].stop();
      if (cf.rank == 0)
        cout << "    Generated in: " << c(cf.colorFlag, CYA) << f->timers["gridTimer"].getf(cf.timeFormat)
             << c(cf.colorFlag, NON) << endl
             << endl;
    }

    if (cf.rank == 0)
      if (cf.write_freq > 0 || cf.restart_freq > 0)
        cout << c(cf.colorFlag, GRE)
             << "Writing Initial Conditions:" << c(cf.colorFlag, NON) << endl;
    if (cf.write_freq > 0) {
      f->timers["solWrite"].reset();
      cf.w->writeSolution(cf, f, 0, 0.00);
      f->timers["solWrite"].accumulate();
    }
    if (cf.restart_freq > 0) {
      f->timers["resWrite"].reset();
      cf.w->writeRestart(cf, f, 0, 0.00);
      f->timers["resWrite"].accumulate();
    }
    f->timers["writeTimer"].stop();
    if (cf.rank == 0)
      if (cf.write_freq > 0 || cf.restart_freq > 0)
        cout << "    Wrote in: " << c(cf.colorFlag, CYA) << f->timers["writeTimer"].getf(cf.timeFormat)
             << c(cf.colorFlag, NON) << endl;

  }else{
    f->timers["writeTimer"].start();
    if (cf.rank == 0)
      cout << c(cf.colorFlag, GRE) << "Loading Restart File:" << c(cf.colorFlag, NON)
           << endl;
    cf.loadTimer.reset();
    cf.w->readSolution(cf, f->grid, f->var);
    cf.loadTimer.stop();
    if (cf.rank == 0)
      cout << "    Loaded in: " << setprecision(2) << c(cf.colorFlag, CYA)
           << cf.loadTimer.getf(cf.timeFormat) << "s" << c(cf.colorFlag, NON) << endl;
    if (cf.rank == 0)
      if (cf.write_freq > 0 || cf.restart_freq > 0)
        cout << c(cf.colorFlag, GRE)
             << "Writing Initial Conditions:" << c(cf.colorFlag, NON) << endl;
    if (cf.write_freq > 0) {
      f->timers["solWrite"].reset();
      cf.w->writeSolution(cf, f, 0, 0.00);
      f->timers["solWrite"].accumulate();
    }
    if (cf.restart_freq > 0) {
      f->timers["resWrite"].reset();
      cf.w->writeRestart(cf, f, 0, 0.00);
      f->timers["resWrite"].accumulate();
    }
    f->timers["writeTimer"].stop();
    if (cf.rank == 0)
      if (cf.write_freq > 0 || cf.restart_freq > 0)
        cout << "    Wrote in: " << c(cf.colorFlag, CYA) << f->timers["writeTimer"].getf(cf.timeFormat)
             << c(cf.colorFlag, NON) << endl;
  }
  // notify simulation start
  if (cf.rank == 0) {
    cout << endl << "-----------------------" << endl << endl;
   cout << c(cf.colorFlag, GRE) << "Starting Simulation:" << c(cf.colorFlag, NON) << endl;
  }
}

void Fiesta::checkIO(struct inputConfig &cf, rk_func *f, int t, double time){
  // Print time step info if necessary
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
           << right << setw(13) << f->timers["gridTimer"].getf(cf.timeFormat) << c(cf.colorFlag, NON) << endl;
      cout << c(cf.colorFlag, NON) << left << setw(36)
           << "    Initial Condition WriteTime:" << c(cf.colorFlag, NON)
           << c(cf.colorFlag, CYA) << right << setw(13) << f->timers["writeTimer"].getf(cf.timeFormat)
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
 
void Fiesta::finalize(){
  Kokkos::finalize();
#ifndef NOMPI
  MPI_Finalize();
#endif
}
