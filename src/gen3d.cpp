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
#include <cassert>
#include "input.hpp"
#ifdef HAVE_MPI
#include "mpi.hpp"
#include "mpi.h"
#endif
#include "Kokkos_Core.hpp"
#include "advect.hpp"
#include "bc.hpp"
#include "debug.hpp"
#include "flux.hpp"
#include "gen3d.hpp"
#include "metric.hpp"
#include "presgrad.hpp"
#include "secondary.hpp"
#include "velocity.hpp"

gen3d_func::gen3d_func(struct inputConfig &cf_) : rk_func(cf_) {

  // Allocate all device arrays here
  grid    = FS4D("coords",  cf.ni,  cf.nj,  cf.nk,  3);      // Grid Coords
  var     = FS4D("var",     cf.ngi, cf.ngj, cf.ngk, cf.nvt); // Primary  Array
  metrics = FS5D("metrics", cf.ngi, cf.ngj, cf.ngk, 3, 3);   // Metric Tensor
  tmp1    = FS4D("tmp1",    cf.ngi, cf.ngj, cf.ngk, cf.nvt); // Temporary Array
  dvar    = FS4D("dvar",    cf.ngi, cf.ngj, cf.ngk, cf.nvt); // RHS Output
  rho     = FS3D("rho",     cf.ngi, cf.ngj, cf.ngk);         // Total Density
  p       = FS3D("p",       cf.ngi, cf.ngj, cf.ngk);         // Pressure
  T       = FS3D("T",       cf.ngi, cf.ngj, cf.ngk);         // Temperature
  tvel    = FS4D("tvel",    cf.ngi, cf.ngj, cf.ngk,  3);     // Velocity
  fluxx   = FS3D("fluxx",   cf.ngi, cf.ngj, cf.ngk); // Advective Fluxes in X
  fluxy   = FS3D("fluxy",   cf.ngi, cf.ngj, cf.ngk); // Advective Fluxes in Y
  fluxz   = FS3D("fluxz",   cf.ngi, cf.ngj, cf.ngk); // Advective Fluxes in z

  // Create and copy minimal configuration array for data needed
  // withing Kokkos kernels.
  cd = FS1D("deviceCF", 6 + cf.ns * 3);
  FS1DH hostcd = Kokkos::create_mirror_view(cd);
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

  // Primaty Variable Names
  varNames.push_back("X-Momentum");
  varNames.push_back("Y-Momentum");
  varNames.push_back("Z-Momentum");
  varNames.push_back("Energy");
  for (int v=0; v<cf.ns; ++v)
    varNames.push_back("Density " + cf.speciesName[v]);
  assert(varNames.size()==cf.nvt);

  // Secondary Variable Names
  varxNames.push_back("X-Velocity");
  varxNames.push_back("Y-Velocity");
  varxNames.push_back("Z-Velocity");
  varxNames.push_back("Pressure");
  varxNames.push_back("Temperature");
  varxNames.push_back("Density");

  // Create Secondary Variable Array
  varx = FS4D("varx",cf.ngi,cf.ngj,cf.ngk,varxNames.size());

  // Create Timers
  timers["flux"]        = Timer::fiestaTimer("Flux Calculation");
  timers["pressgrad"]   = Timer::fiestaTimer("Pressure Gradient Calculation");
  timers["calcSecond"]  = Timer::fiestaTimer("Secondary Variable Calculation");
  timers["solWrite"]    = Timer::fiestaTimer("Solution Write Time");
  timers["resWrite"]    = Timer::fiestaTimer("Restart Write Time");
  timers["statCheck"]   = Timer::fiestaTimer("Status Check");
  timers["rk"] = Timer::fiestaTimer("Runge Stage Update");
  timers["halo"] = Timer::fiestaTimer("Halo Exchanges");
  timers["bc"] = Timer::fiestaTimer("Boundary Conditions");
  timers["calcMetrics"] = Timer::fiestaTimer("Metric Computations");
};

void gen3d_func::preStep() {}

void gen3d_func::postStep() {
  
  // MDRange Policy for all cells including ghost cells
  policy_f3 ghost_pol = policy_f3({0, 0, 0}, {cf.ngi, cf.ngj, cf.ngk});

  // Copy secondary variables to extra variables array
  if (( (cf.write_freq >0) && (cf.t % cf.write_freq == 0) )||
      ( (cf.stat_freq  >0) && (cf.t % cf.stat_freq  == 0) )){

      timers["calcSecond"].reset();
      Kokkos::parallel_for(ghost_pol, calculateRhoPT3D(var, p, rho, T, cd));
      Kokkos::parallel_for(ghost_pol, computeVelocity3D(var, rho, tvel));
      Kokkos::parallel_for(ghost_pol, copyExtraVars3D(varx, tvel, p, rho, T));
      Kokkos::fence();
      timers["calcSecond"].accumulate();

    }
}
void gen3d_func::preSim() {
  policy_f3 cell_pol3 = policy_f3(
      {cf.ng, cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng, cf.ngk - cf.ng});

  // compute metrics
  timers["calcMetrics"].reset();
  Kokkos::parallel_for(cell_pol3, computeMetrics3D(metrics, grid));

  FS5D &mt = metrics;

  int mngi = cf.ngi;
  int mngj = cf.ngj;
  int mngk = cf.ngk;
  int mnci = cf.nci;
  int mncj = cf.ncj;
  int mnck = cf.nck;
  int mng = cf.ng;

#ifdef HAVE_MPI
  // mpi exchange of metrics
  ls = FS5D("leftSend", cf.ng, cf.ngj, cf.ngk, 3, 3);
  lr = FS5D("leftRecv", cf.ng, cf.ngj, cf.ngk, 3, 3);
  rs = FS5D("rightSend", cf.ng, cf.ngj, cf.ngk, 3, 3);
  rr = FS5D("rightRecv", cf.ng, cf.ngj, cf.ngk, 3, 3);
  bs = FS5D("bottomSend", cf.ngi, cf.ng, cf.ngk, 3, 3);
  br = FS5D("bottomRecv", cf.ngi, cf.ng, cf.ngk, 3, 3);
  ts = FS5D("topSend", cf.ngi, cf.ng, cf.ngk, 3, 3);
  tr = FS5D("topRecv", cf.ngi, cf.ng, cf.ngk, 3, 3);
  hs = FS5D("hindSend", cf.ngi, cf.ngj, cf.ng, 3, 3);
  hr = FS5D("hindRecv", cf.ngi, cf.ngj, cf.ng, 3, 3);
  fs = FS5D("frontSend", cf.ngi, cf.ngj, cf.ng, 3, 3);
  fr = FS5D("frontRecv", cf.ngi, cf.ngj, cf.ng, 3, 3);

  lsH = Kokkos::create_mirror_view(ls);
  lrH = Kokkos::create_mirror_view(lr);
  rsH = Kokkos::create_mirror_view(rs);
  rrH = Kokkos::create_mirror_view(rr);
  bsH = Kokkos::create_mirror_view(bs);
  brH = Kokkos::create_mirror_view(br);
  tsH = Kokkos::create_mirror_view(ts);
  trH = Kokkos::create_mirror_view(tr);
  hsH = Kokkos::create_mirror_view(hs);
  hrH = Kokkos::create_mirror_view(hr);
  fsH = Kokkos::create_mirror_view(fs);
  frH = Kokkos::create_mirror_view(fr);

  policy_f5 xPol = policy_f5({0, 0, 0, 0, 0}, {cf.ng, cf.ngj, cf.ngk, 3, 3});
  policy_f5 yPol = policy_f5({0, 0, 0, 0, 0}, {cf.ngi, cf.ng, cf.ngk, 3, 3});
  policy_f5 zPol = policy_f5({0, 0, 0, 0, 0}, {cf.ngi, cf.ngj, cf.ng, 3, 3});

  FS5D &lsR = ls;
  FS5D &rsR = rs;
  FS5D &bsR = bs;
  FS5D &tsR = ts;
  FS5D &hsR = hs;
  FS5D &fsR = fs;
  Kokkos::parallel_for(
      xPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int m,
                          const int n) {
        lsR(i, j, k, m, n) = mt(mng + i, j, k, m, n);
        rsR(i, j, k, m, n) = mt(mnci + i, j, k, m, n);
      });
  Kokkos::parallel_for(
      yPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int m,
                          const int n) {
        bsR(i, j, k, m, n) = mt(i, mng + j, k, m, n);
        tsR(i, j, k, m, n) = mt(i, mncj + j, k, m, n);
      });
  Kokkos::parallel_for(
      zPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int m,
                          const int n) {
        hsR(i, j, k, m, n) = mt(i, j, mng + k, m, n);
        fsR(i, j, k, m, n) = mt(i, j, mnck + k, m, n);
      });
  Kokkos::deep_copy(lsH, ls);
  Kokkos::deep_copy(rsH, rs);
  Kokkos::deep_copy(bsH, bs);
  Kokkos::deep_copy(tsH, ts);
  Kokkos::deep_copy(hsH, hs);
  Kokkos::deep_copy(fsH, fs);

  MPI_Request reqs[12];

  MPI_Isend(lsH.data(), cf.ng * cf.ngj * cf.ngk * 9, MPI_DOUBLE, cf.xMinus, 0,
            cf.comm, &reqs[0]);
  MPI_Irecv(lrH.data(), cf.ng * cf.ngj * cf.ngk * 9, MPI_DOUBLE, cf.xMinus,
            MPI_ANY_TAG, cf.comm, &reqs[1]);

  MPI_Isend(rsH.data(), cf.ng * cf.ngj * cf.ngk * 9, MPI_DOUBLE, cf.xPlus, 0,
            cf.comm, &reqs[2]);
  MPI_Irecv(rrH.data(), cf.ng * cf.ngj * cf.ngk * 9, MPI_DOUBLE, cf.xPlus,
            MPI_ANY_TAG, cf.comm, &reqs[3]);

  MPI_Isend(bsH.data(), cf.ngi * cf.ng * cf.ngk * 9, MPI_DOUBLE, cf.yMinus, 0,
            cf.comm, &reqs[4]);
  MPI_Irecv(brH.data(), cf.ngi * cf.ng * cf.ngk * 9, MPI_DOUBLE, cf.yMinus,
            MPI_ANY_TAG, cf.comm, &reqs[5]);

  MPI_Isend(tsH.data(), cf.ngi * cf.ng * cf.ngk * 9, MPI_DOUBLE, cf.yPlus, 0,
            cf.comm, &reqs[6]);
  MPI_Irecv(trH.data(), cf.ngi * cf.ng * cf.ngk * 9, MPI_DOUBLE, cf.yPlus,
            MPI_ANY_TAG, cf.comm, &reqs[7]);

  MPI_Isend(hsH.data(), cf.ngi * cf.ngj * cf.ng * 9, MPI_DOUBLE, cf.zMinus, 0,
            cf.comm, &reqs[8]);
  MPI_Irecv(hrH.data(), cf.ngi * cf.ngj * cf.ng * 9, MPI_DOUBLE, cf.zMinus,
            MPI_ANY_TAG, cf.comm, &reqs[9]);

  MPI_Isend(fsH.data(), cf.ngi * cf.ngj * cf.ng * 9, MPI_DOUBLE, cf.zPlus, 0,
            cf.comm, &reqs[10]);
  MPI_Irecv(frH.data(), cf.ngi * cf.ngj * cf.ng * 9, MPI_DOUBLE, cf.zPlus,
            MPI_ANY_TAG, cf.comm, &reqs[11]);

  MPI_Waitall(12, reqs, MPI_STATUS_IGNORE);

  Kokkos::deep_copy(lr, lrH);
  Kokkos::deep_copy(rr, rrH);
  Kokkos::deep_copy(br, brH);
  Kokkos::deep_copy(tr, trH);
  Kokkos::deep_copy(hr, hrH);
  Kokkos::deep_copy(fr, frH);
  FS5D &lrR = lr;
  FS5D &rrR = rr;
  FS5D &brR = br;
  FS5D &trR = tr;
  FS5D &hrR = hr;
  FS5D &frR = fr;

  Kokkos::parallel_for(
      xPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int m,
                          const int n) {
        mt(i, j, k, m, n) = lrR(i, j, k, m, n);
        mt(mngi - mng + i, j, k, m, n) = rrR(i, j, k, m, n);
      });
  Kokkos::parallel_for(
      yPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int m,
                          const int n) {
        mt(i, j, k, m, n) = brR(i, j, k, m, n);
        mt(i, mngj - mng + j, k, m, n) = trR(i, j, k, m, n);
      });
  Kokkos::parallel_for(
      zPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int m,
                          const int n) {
        mt(i, j, k, m, n) = hrR(i, j, k, m, n);
        mt(i, j, mngk - mng + k, m, n) = frR(i, j, k, m, n);
      });

#endif

  if (cf.xMinus < 0) {
    Kokkos::parallel_for(
        "metricbcl", policy_f3({0, 0, 0}, {mng, mngj, mngk}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
              mt(mng - i - 1, j, k, m, n) = mt(mng + i, j, k, m, n);
        });
  }
  if (cf.xPlus < 0) {
    Kokkos::parallel_for(
        "metricbcr", policy_f3({0, 0, 0}, {mng, mngj, mngk}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
              mt(mngi - 1 - i, j, k, m, n) = mt(mnci + i, j, k, m, n);
        });
  }
  if (cf.yMinus < 0) {
    Kokkos::parallel_for(
        "metricbcb", policy_f3({0, 0, 0}, {mngi, mng, mngk}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
              mt(i, mng - j - 1, k, m, n) = mt(i, mng + j, k, m, n);
        });
  }
  if (cf.yPlus < 0) {
    Kokkos::parallel_for(
        "metricbct", policy_f3({0, 0, 0}, {mngi, mng, mngk}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
              mt(i, mngj - 1 - j, k, m, n) = mt(i, mncj + j, k, m, n);
        });
  }
  if (cf.zMinus < 0) {
    Kokkos::parallel_for(
        "metricbch", policy_f3({0, 0, 0}, {mngi, mngj, mng}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
              mt(i, j, mng - k - 1, m, n) = mt(i, j, mng + k, m, n);
        });
  }
  if (cf.zPlus < 0) {
    Kokkos::parallel_for(
        "metricbcf", policy_f3({0, 0, 0}, {mngi, mngj, mng}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
              mt(i, j, mngk - 1 - k, m, n) = mt(i, j, mnck + k, m, n);
        });
  }

  Kokkos::fence();
  timers["calcMetrics"].accumulate();
}
void gen3d_func::postSim() {}

void gen3d_func::compute() {
  // create range policies
  policy_f3 ghost_pol = policy_f3({0, 0, 0}, {cf.ngi, cf.ngj, cf.ngk});
  policy_f3 cell_pol = policy_f3(
      {cf.ng, cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng, cf.ngk - cf.ng});
  policy_f3 weno_pol =
      policy_f3({cf.ng - 1, cf.ng - 1, cf.ng - 1},
                {cf.ngi - cf.ng, cf.ngj - cf.ng, cf.ngk - cf.ng});

  /**** WENO ****/
  timers["calcSecond"].reset();
  Kokkos::parallel_for(ghost_pol, calculateRhoPT3D(var, p, rho, T, cd));
  Kokkos::parallel_for(ghost_pol,
                       computeGenVelocity3D(var, metrics, rho, tvel));
  Kokkos::fence();
  timers["calcSecond"].accumulate();

  timers["flux"].reset();
  for (int v = 0; v < cf.nv; ++v) {
    Kokkos::parallel_for(
        weno_pol,
        calculateFluxesG(var, p, rho, tvel, fluxx, fluxy, fluxz, cf.dx, cf.dy, cf.dx, v));
    Kokkos::parallel_for(cell_pol, advect3D(dvar, varx, fluxx, fluxy, fluxz, v));
  }
  Kokkos::fence();
  timers["flux"].accumulate();

  timers["pressgrad"].reset();
  Kokkos::parallel_for(cell_pol,
                       applyGenPressureGradient3D(dvar, metrics, p, cd));
  Kokkos::fence();
  timers["pressgrad"].accumulate();
}
