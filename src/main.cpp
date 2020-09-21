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

using namespace std;

int main(int argc, char *argv[]) {
  

  struct inputConfig cf = Fiesta::initialize(argc,argv);

  cf.totalTimer.start();
  cf.initTimer.start();

  {
#ifndef NOMPI
    // create mpi buffers
    mpiBuffers m(cf);
#endif

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

    Fiesta::initializeSimulation(cf,f);

    // execution policy for Runge Kutta update kernels, should be moved to
    // header
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_1;

    // Variables to track time step information
    double time = cf.time;
    int tstart = cf.tstart;

    // pre simulation hook
    f->preSim();

    // stop initialization timer
    cf.initTimer.stop();

    // reset simulation timer
    cf.simTimer.reset();

    // notify simulation start
    if (cf.rank == 0) {
      cout << endl << "-----------------------" << endl << endl;
      cout << c(cf.colorFlag, GRE) << "Starting Simulation:" << c(cf.colorFlag, NON) << endl;
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
      applyBCs(cf, f->var);

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
      applyBCs(cf, f->var);

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

      Fiesta::checkIO(cf,f,t,time);
    } // End main time loop

    // post simulation hook
    f->postSim();

    // stop simulation timer
    cf.simTimer.stop();

    // notify simulation complete
    if (cf.rank == 0)
      cout << c(cf.colorFlag, GRE) << "Simulation Complete!" << c(cf.colorFlag, NON) << endl;

    Fiesta::reportTimers(cf,f);
  }

  Fiesta::finalize();

  return 0;
}
