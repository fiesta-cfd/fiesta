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

using namespace std;

// void fnExit1(void){
//    Kokkos::finalize();
//}

int main(int argc, char *argv[]) {
  printf("LUA TEST\n");
  luaReader L("fiesta.lua");

  double nd = L.getDouble("dx");
  printf("dx: %f\n",nd);

  double dd = L.getDouble("dr",0.123);
  printf("dr: %f\n",dd);

  int    ni = L.getInt("ni");
  printf("ni: %d\n",ni);

  int    di = L.getInt("nm",457);
  printf("nm: %d\n",di);

  std::string ns = L.getString("title");
  printf("title: %s\n",ns.c_str());

  std::string ds = L.getString("kokkos","generalized");
  printf("kokkos: %s\n",ds.c_str());

  bool nb = L.getBool("visc");
  printf("visc: %d\n",nb);

  bool db = L.getBool("fiesta",false);
  printf("fiesta: %d\n",db);

  double *Ma = (double*)malloc(2*sizeof(double));
  L.getDoubles("mu",2,Ma);
  printf("mu: %f, %f\n",Ma[0],Ma[1]);

  int *Mb = (int*)malloc(2*sizeof(int));
  L.getInts("luna",2,Mb);
  printf("luna: %d, %d\n",Mb[0],Mb[1]);

  std::vector<std::string> Mc;
  printf("#######\n");
  L.getStrings("sol",3,Mc);
  printf("#######\n");
  printf("sol: %s %s %s\n",Mc[0].c_str(),Mc[1].c_str(),Mc[2].c_str());

  L.close();
  printf("LUA TEST\n");
  

  // Get the command line options including input file name if supplied.
  struct commandArgs cArgs = getCommandlineOptions(argc, argv);

  // If the --version command line option is present, then just print version
  // and compilation information and exit.
  if (cArgs.versionFlag) {
    cout << "Fiesta" << endl;
    cout << "Version:    '" << FIESTA_VERSION << "'" << endl;
    cout << "Build Type: '" << FIESTA_OPTIONS << "'" << endl;
    cout << "Build Time: '" << FIESTA_BTIME << "'" << endl;
    exit(EXIT_SUCCESS);
  }

  // Create and start some timers
  fiestaTimer totalTimer;    // Total program timer
  fiestaTimer initTimer;     // Timer for program initialization
  fiestaTimer solWriteTimer; // Timer for solution file write times
  fiestaTimer simTimer;      // Timer for simulation time steppin
  fiestaTimer loadTimer;     // Timer for I.C. generation or restart read
  fiestaTimer resWriteTimer; // Timer for restart file writes
  fiestaTimer gridTimer;     // Timer for grid generation
  fiestaTimer writeTimer;    // Timer for

  // Initialize command MPI and get rank.
  int temp_rank;
  temp_rank = 0;
#ifndef NOMPI
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &temp_rank);
#endif

  // Print spash screen with logo, version number and compilation info
  // to signify process startup.
  if (temp_rank == 0)
    printSplash(cArgs.colorFlag);

  // Initialize kokkos and set kokkos finalize as exit function.
  Kokkos::initialize(argc, argv);
  //    atexit(fnExit1);

  // Execute input file and generate simulation configuration
  struct inputConfig cf;
  cf = executeConfiguration(cArgs);

  {
#ifndef NOMPI
    // Perform MPI decomposition and assign cells to ranks
    cf = mpi_init(cf);

    // create mpi buffers
    mpiBuffers m(cf);
#endif

    // Print out configuration
    if (cf.rank == 0)
      printConfig(cf, cArgs.colorFlag);

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

    /*** Choose Scheme ***/
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

    // Create writer object
#ifdef NOMPI
    serialVTKWriter w(cf, f->grid, f->var);
#else
    fstWriter w(cf, f);
    // cgnsWriter w(cf,f->grid,f->var);
#endif

    // If not restarting, generate initial conditions and grid
    if (cf.restart == 0) {
      if (cf.rank == 0)
        cout << c(cArgs.colorFlag, GRE)
             << "Generating Initial Conditions:" << c(cArgs.colorFlag, NON) << endl;
      loadTimer.start();
      loadInitialConditions(cf, f->var);
      loadTimer.stop();
      if (cf.rank == 0)
        cout << "    Generated in: " << c(cArgs.colorFlag, CYA) << loadTimer.getf(cArgs.timeFormat)
             << c(cArgs.colorFlag, NON) << endl
             << endl;

      if (cf.rank == 0)
        cout << c(cArgs.colorFlag, GRE) << "Generating Grid:" << c(cArgs.colorFlag, NON) << endl;
      gridTimer.start();
      loadGrid(cf, f->grid);
      gridTimer.stop();
      if (cf.rank == 0)
        cout << "    Generated in: " << c(cArgs.colorFlag, CYA) << gridTimer.getf(cArgs.timeFormat)
             << c(cArgs.colorFlag, NON) << endl
             << endl;

      if (cf.particle == 1) {
        if (cf.rank == 0)
          cout << c(cArgs.colorFlag, GRE) << "Generating Particles:" << c(cArgs.colorFlag, NON)
               << endl;
        gridTimer.start();
        loadParticles(cf, f->particles);
        gridTimer.stop();
        if (cf.rank == 0)
          cout << "    Generated in: " << c(cArgs.colorFlag, CYA) << gridTimer.getf(cArgs.timeFormat)
               << c(cArgs.colorFlag, NON) << endl
               << endl;
      }
    }

    // execution policy for Runge Kutta update kernels, should be moved to
    // header
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_1;

    // Variables to track time step information
    double time = cf.time;
    int tstart = cf.tstart;

    // Read Restart or Write initial conditions
    //#ifndef NOMPI
    writeTimer.start();
    if (cf.restart == 1) {
      if (cf.rank == 0)
        cout << c(cArgs.colorFlag, GRE) << "Loading Restart File:" << c(cArgs.colorFlag, NON)
             << endl;
      loadTimer.reset();
      w.readSolution(cf, f->grid, f->var);
      loadTimer.stop();
      if (cf.rank == 0)
        cout << "    Loaded in: " << setprecision(2) << c(cArgs.colorFlag, CYA)
             << loadTimer.getf(cArgs.timeFormat) << "s" << c(cArgs.colorFlag, NON) << endl;
      // if (cf.rank == 0)
      //   if (cf.write_freq > 0 || cf.restart_freq > 0)
      //     cout << c(cArgs.colorFlag, GRE)
      //          << "Writing Initial Conditions:" << c(cArgs.colorFlag, NON) << endl;
      // if (cf.write_freq > 0) {
      //   solWriteTimer.reset();
      //   w.writeSolution(cf, f, 0, 0.00);
      //   solWriteTimer.accumulate();
      // }
      // writeTimer.stop();
      // if (cf.rank == 0)
      //   if (cf.write_freq > 0 || cf.restart_freq > 0)
      //     cout << "    Wrote in: " << c(cArgs.colorFlag, CYA) << writeTimer.getf(cArgs.timeFormat)
      //          << c(cArgs.colorFlag, NON) << endl;
    } else {
      if (cf.rank == 0)
        if (cf.write_freq > 0 || cf.restart_freq > 0)
          cout << c(cArgs.colorFlag, GRE)
               << "Writing Initial Conditions:" << c(cArgs.colorFlag, NON) << endl;
      if (cf.write_freq > 0) {
        solWriteTimer.reset();
        w.writeSolution(cf, f, 0, 0.00);
        solWriteTimer.accumulate();
      }
      if (cf.restart_freq > 0) {
        resWriteTimer.reset();
        w.writeRestart(cf, f, 0, 0.00);
        resWriteTimer.accumulate();
      }
      writeTimer.stop();
      if (cf.rank == 0)
        if (cf.write_freq > 0 || cf.restart_freq > 0)
          cout << "    Wrote in: " << c(cArgs.colorFlag, CYA) << writeTimer.getf(cArgs.timeFormat)
               << c(cArgs.colorFlag, NON) << endl;
    }

    //#endif

    // pre simulation hook
    f->preSim();

    // stop initialization timer
    initTimer.stop();

    // reset simulation timer
    simTimer.reset();

    // notify simulation start
    if (cf.rank == 0) {
      cout << endl << "-----------------------" << endl << endl;
      cout << c(cArgs.colorFlag, GRE) << "Starting Simulation:" << c(cArgs.colorFlag, NON) << endl;
    }

    // // // // // // // //  \\ \\ \\ \\ \\ \\ \\ \\
 // // // // // // MAIN TIME LOOP \\ \\ \\ \\ \\ \\
// // // // // // // // //\\ \\ \\ \\ \\ \\ \\ \\ \\

    for (int t = tstart; t < cf.tend; ++t) {
      time = time + cf.dt;
      cf.t = t + 1;

      // pre timestep hook
      f->preStep();

      // apply boundary conditions
#ifndef NOMPI
      applyBCs(cf, f->var, m);
#else
      applyBCs(cf, f->var);
#endif

      // First Stage Compute
      f->compute();

      // assign temporary variables
      FS4D mytmp = f->tmp1;
      FS4D myvar = f->var;
      FS4D mydvar = f->dvar;

      // First stage update
      f->timers["rk"].reset();
      Kokkos::parallel_for(
          "Loop1", policy_1({0, 0, 0}, {cf.ngi, cf.ngj, cf.ngk}),
          KOKKOS_LAMBDA(const int i, const int j, const int k) {
            for (int v=0; v<cf.nvt; ++v){
              mytmp(i, j, k, v) = myvar(i, j, k, v);
              myvar(i, j, k, v) =
                myvar(i, j, k, v) + 0.5 * cf.dt * mydvar(i, j, k, v);
            }
          });
      Kokkos::fence();
      f->timers["rk"].accumulate();

      // apply boundary conditions
#ifndef NOMPI
      applyBCs(cf, f->var, m);
#else
      applyBCs(cf, f->var);
#endif

      // Second stage compute
      f->compute();

      // assign temporary variables
      // mytmp = f->tmp1;
      myvar = f->var;
      mydvar = f->dvar;

      // Second stage update
      f->timers["rk"].reset();
      Kokkos::parallel_for(
          "Loop2", policy_1({0, 0, 0}, {cf.ngi, cf.ngj, cf.ngk}),
          KOKKOS_LAMBDA(const int i, const int j, const int k) {
            for (int v=0; v<cf.nvt; ++v){
              myvar(i, j, k, v) = mytmp(i, j, k, v) + cf.dt * mydvar(i, j, k, v);
            }
          });
      Kokkos::fence();
      f->timers["rk"].accumulate();

      // post timestep hook
      f->postStep();

      // Print time step info if necessary
      if (cf.rank == 0) {
        if (cf.out_freq > 0)
          if ((t + 1) % cf.out_freq == 0)
            cout << c(cArgs.colorFlag, YEL) << left << setw(15)
                 << "    Iteration:" << c(cArgs.colorFlag, NON) << c(cArgs.colorFlag, CYA) << right
                 << setw(0) << t + 1 << c(cArgs.colorFlag, NON) << "/" << c(cArgs.colorFlag, CYA)
                 << left << setw(0) << cf.tend << c(cArgs.colorFlag, NON) << ", "
                 << c(cArgs.colorFlag, CYA) << right << setw(0) << setprecision(3)
                 << scientific << time << "s" << c(cArgs.colorFlag, NON) << endl;
      }
      //#ifndef NOMPI
      // Write solution file if necessary
      if (cf.write_freq > 0) {
        if ((t + 1) % cf.write_freq == 0) {
          f->timers["solWrite"].reset();
          w.writeSolution(cf, f, t + 1, time);
          Kokkos::fence();
          f->timers["solWrite"].accumulate();
        }
      }
      // Write restart file if necessary
      if (cf.restart_freq > 0) {
        if ((t + 1) % cf.restart_freq == 0) {
          f->timers["resWrite"].reset();
          w.writeRestart(cf, f, t + 1, time);
          Kokkos::fence();
          f->timers["resWrite"].accumulate();
        }
      }
      //#endif
      // Print status check if necessary
      if (cf.stat_freq > 0) {
        if ((t + 1) % cf.stat_freq == 0) {
          f->timers["statCheck"].reset();
          statusCheck(cArgs.colorFlag, cf, f, time, totalTimer, simTimer);
          f->timers["statCheck"].accumulate();
        }
      }
    } // End main time loop

    // post simulation hook
    f->postSim();

    // stop simulation timer
    simTimer.stop();

    // notify simulation complete
    if (cf.rank == 0)
      cout << c(cArgs.colorFlag, GRE) << "Simulation Complete!" << c(cArgs.colorFlag, NON) << endl;

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

    // stop total program timer
    totalTimer.stop();

    // prin timer values
    if (cf.rank == 0) {
      cout << endl << "-----------------------" << endl << endl;
      cout.precision(2);
      cout << c(cArgs.colorFlag, GRE) << left << setw(36)
           << "Total Time:" << c(cArgs.colorFlag, NON) << c(cArgs.colorFlag, CYA) << right
           << setw(13) << totalTimer.getf(cArgs.timeFormat) << c(cArgs.colorFlag, NON) << endl
           << endl;

      cout << c(cArgs.colorFlag, GRE) << left << setw(36)
           << "  Setup Time:" << c(cArgs.colorFlag, NON) << c(cArgs.colorFlag, CYA) << right
           << setw(13) << initTimer.getf(cArgs.timeFormat) << c(cArgs.colorFlag, NON) << endl;
      if (cf.restart == 1)
        cout << c(cArgs.colorFlag, NON) << left << setw(36)
             << "    Restart Read:" << c(cArgs.colorFlag, NON) << c(cArgs.colorFlag, CYA) << right
             << setw(13) << loadTimer.getf(cArgs.timeFormat) << c(cArgs.colorFlag, NON) << endl
             << endl;
      else {
        cout << c(cArgs.colorFlag, NON) << left << setw(36)
             << "    Initial Condition Generation:" << c(cArgs.colorFlag, NON)
             << c(cArgs.colorFlag, CYA) << right << setw(13) << loadTimer.getf(cArgs.timeFormat)
             << c(cArgs.colorFlag, NON) << endl;
        cout << c(cArgs.colorFlag, NON) << left << setw(36)
             << "    Grid Generation:" << c(cArgs.colorFlag, NON) << c(cArgs.colorFlag, CYA)
             << right << setw(13) << gridTimer.getf(cArgs.timeFormat) << c(cArgs.colorFlag, NON) << endl;
        cout << c(cArgs.colorFlag, NON) << left << setw(36)
             << "    Initial Consition WriteTime:" << c(cArgs.colorFlag, NON)
             << c(cArgs.colorFlag, CYA) << right << setw(13) << writeTimer.getf(cArgs.timeFormat)
             << c(cArgs.colorFlag, NON) << endl
             << endl;
      }

      cout << c(cArgs.colorFlag, GRE) << left << setw(36)
           << "  Simulation Time:" << c(cArgs.colorFlag, NON) << c(cArgs.colorFlag, CYA) << right
           << setw(13) << simTimer.getf(cArgs.timeFormat) << c(cArgs.colorFlag, NON) << endl;
      for (auto tmr : stmr) {
        cout << c(cArgs.colorFlag, NON) << left << setw(36)
             << "    " + tmr.second.describe() + ":" << c(cArgs.colorFlag, NON)
             << c(cArgs.colorFlag, CYA) << right << setw(13) << tmr.second.getf(cArgs.timeFormat)
             << c(cArgs.colorFlag, NON) << endl;
      }
      cout << " " << endl;
    }

    // delete &m;
  }

  Kokkos::finalize();
  // finalize mpi
#ifndef NOMPI
  MPI_Finalize();
#endif

  // exit success
  return 0;
}
