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

#include "gen2d.hpp"
#include <cassert>

gen2d_func::gen2d_func(struct inputConfig &cf_) : rk_func(cf_) {
  // Allocate all device variables here
  grid    = FS4D("grid",    cf.ni,  cf.nj,  cf.nk, 3);       // Grid Coords
  metrics = FS4D("metrics", cf.ngi, cf.ngj, 2, 2);           // Metric Tensor
  var     = FS4D("var",     cf.ngi, cf.ngj, cf.ngk, cf.nvt); // Primary Vars
  tmp1    = FS4D("tmp1",    cf.ngi, cf.ngj, cf.ngk, cf.nvt); // Temporary Vars
  dvar    = FS4D("dvar",    cf.ngi, cf.ngj, cf.ngk, cf.nvt); // RHS Output
  tvel    = FS3D("vel",     cf.ngi, cf.ngj, 2);              // Velocity
  p       = FS2D("p",       cf.ngi, cf.ngj);                 // Pressure
  T       = FS2D("T",       cf.ngi, cf.ngj);                 // Temperature
  rho     = FS2D("rho",     cf.ngi, cf.ngj);                 // Total Density
  fluxx   = FS2D("fluxx",   cf.ngi, cf.ngj);              // Advective Fluxes X
  fluxy   = FS2D("fluxy",   cf.ngi, cf.ngj);              // Advective Fluxes Y

  // Primary Variable Names
  varNames.push_back("X-Momentum");
  varNames.push_back("Y-Momentum");
  varNames.push_back("Energy");
  for (int v=0; v<cf.ns; ++v)
    varNames.push_back("Density " + cf.speciesName[v]);
  assert(varNames.size() == cf.nvt);

  if (cf.noise == 1) {
    noise = FS2D_I("noise", cf.ngi, cf.ngj); // Noise indicator array
  }

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

  // Secondary Variable Names
  varxNames.push_back("X-Velocity");
  varxNames.push_back("Y-Velocity");
  varxNames.push_back("Pressure");
  varxNames.push_back("Temperature");
  varxNames.push_back("Density");

  varx  = FS4D("varx", cf.ngi, cf.ngj, cf.ngk, varxNames.size()); // Extra Vars

  // Create Timers
  timers["flux"] = Timer::fiestaTimer("Flux Calculation");
  timers["pressgrad"] = Timer::fiestaTimer("Pressure Gradient Calculation");
  timers["calcSecond"] = Timer::fiestaTimer("Secondary Variable Calculation");
  timers["solWrite"] = Timer::fiestaTimer("Solution Write Time");
  timers["resWrite"] = Timer::fiestaTimer("Restart Write Time");
  timers["statCheck"] = Timer::fiestaTimer("Status Check");
  timers["rk"] = Timer::fiestaTimer("Runge Stage Update");
  timers["halo"] = Timer::fiestaTimer("Halo Exchanges");
  timers["bc"] = Timer::fiestaTimer("Boundary Conditions");
  timers["calcMetrics"] = Timer::fiestaTimer("Metric Calculations");
  if (cf.noise == 1) {
    timers["noise"] = Timer::fiestaTimer("Noise Removal");
  }

  ghostPol = policy_f({0, 0}, {cf.ngi, cf.ngj});
  cellPol = policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});
  facePol = policy_f({cf.ng - 1, cf.ng - 1}, {cf.ngi - cf.ng, cf.ngj - cf.ng});
};

void gen2d_func::compute() {
  using namespace Kokkos;

  // Calcualte Total Density and Pressure Fields
  timers["calcSecond"].reset();
  parallel_for(ghostPol, calculateRhoPT2D(var, p, rho, T, cd));
  parallel_for(ghostPol, computeGenVelocity2D(var, metrics, rho, tvel));
  fence();
  timers["calcSecond"].accumulate();

  // Calculate and apply weno fluxes for each variable
  timers["flux"].reset();
  for (int v = 0; v < cf.nv; ++v) {
    if (cf.scheme == 3) {
      parallel_for(facePol, computeFluxQuick2D(var,p,rho,fluxx,fluxy,cd,v));
    } else if (cf.scheme == 2) {
      parallel_for(facePol, computeFluxCentered2D(var,p,rho,fluxx,fluxy,cd,v));
    } else {
      parallel_for(facePol, computeFluxWeno2D(var,p,rho,tvel,fluxx,fluxy,cf.dx,cf.dy,v));
    }
    parallel_for(cellPol, advect2D(dvar, fluxx, fluxy, cd, v));
  }
  fence();
  timers["flux"].accumulate();

  // Apply Pressure Gradient Term
  timers["pressgrad"].reset();
  parallel_for(cellPol, applyGenPressureGradient2D(dvar, metrics, p, cd));
  fence();
  timers["pressgrad"].accumulate();
}

void gen2d_func::preStep() {
  /* Pre time step hook executed before the first runge kutta stage */
}

void gen2d_func::postStep() {

  if (( (cf.write_freq >0) && (cf.t % cf.write_freq == 0) )||
      ( (cf.stat_freq  >0) && (cf.t % cf.stat_freq  == 0) )){

      timers["calcSecond"].reset();
      Kokkos::parallel_for(ghostPol, calculateRhoPT2D(var, p, rho, T, cd));
      Kokkos::parallel_for(ghostPol, computeVelocity2D(var, rho, tvel));
      Kokkos::parallel_for(ghostPol, copyExtraVars2D(varx, tvel, p, rho, T));
      Kokkos::fence();
      timers["calcSecond"].accumulate();
    }
  

  if (cf.noise == 1) {
    int M = 0;
    int N = 0;
    FSCAL coff;

    if ((cf.nci - 1) % 2 == 0)
      M = (cf.nci - 1) / 2;
    else
      M = cf.nci / 2;

    if ((cf.ncj - 1) % 2 == 0)
      N = (cf.ncj - 1) / 2;
    else
      N = cf.ncj / 2;
    policy_f noise_pol = policy_f({0, 0}, {M, N});

    coff = 0.0;

    timers["noise"].reset();
    for (int v = 0; v < 2; ++v) {
      Kokkos::parallel_for(noise_pol,
                           detectNoise2D(var, varx, noise, cf.n_dh, coff, cd, v));
      for (int tau = 0; tau < cf.n_nt; ++tau) {
        Kokkos::parallel_for(cellPol, removeNoise2D(dvar, var, varx, noise, cf.n_eta, cd, v));
        //Kokkos::parallel_for(cellPol, updateNoise2D(dvar, var, v));
      }
    }
    Kokkos::fence();
    timers["noise"].accumulate();
  } // end noise
}

void gen2d_func::preSim() {

  // execution policy for all cells including ghost cells
  policy_f cell_pol =
      policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});

  // compute metrics
  timers["calcMetrics"].reset();
  Kokkos::parallel_for(cell_pol, computeMetrics2D(metrics, grid));

  FS4D &mt = metrics;

  int mng = cf.ng;
  int mngi = cf.ngi;
  int mngj = cf.ngj;
  int mnci = cf.nci;
  int mncj = cf.ncj;

#ifdef HAVE_MPI
  // mpi exchange of metrics
  ls = Kokkos::View<FSCAL ****, FS_LAYOUT>("leftSend", cf.ng, cf.ngj, 2, 2);
  lr = Kokkos::View<FSCAL ****, FS_LAYOUT>("leftRecv", cf.ng, cf.ngj, 2, 2);
  rs = Kokkos::View<FSCAL ****, FS_LAYOUT>("rightSend", cf.ng, cf.ngj, 2, 2);
  rr = Kokkos::View<FSCAL ****, FS_LAYOUT>("rightRecv", cf.ng, cf.ngj, 2, 2);
  bs = Kokkos::View<FSCAL ****, FS_LAYOUT>("bottomSend", cf.ngi, cf.ng, 2, 2);
  br = Kokkos::View<FSCAL ****, FS_LAYOUT>("bottomRecv", cf.ngi, cf.ng, 2, 2);
  ts = Kokkos::View<FSCAL ****, FS_LAYOUT>("topSend", cf.ngi, cf.ng, 2, 2);
  tr = Kokkos::View<FSCAL ****, FS_LAYOUT>("topRecv", cf.ngi, cf.ng, 2, 2);

  lsH = Kokkos::create_mirror_view(ls);
  lrH = Kokkos::create_mirror_view(lr);
  rsH = Kokkos::create_mirror_view(rs);
  rrH = Kokkos::create_mirror_view(rr);
  bsH = Kokkos::create_mirror_view(bs);
  brH = Kokkos::create_mirror_view(br);
  tsH = Kokkos::create_mirror_view(ts);
  trH = Kokkos::create_mirror_view(tr);

  FS4D &lsR = ls;
  FS4D &lrR = lr;
  FS4D &rsR = rs;
  FS4D &rrR = rr;
  FS4D &bsR = bs;
  FS4D &brR = br;
  FS4D &tsR = ts;
  FS4D &trR = tr;

  policy_f4 xPol = policy_f4({0, 0, 0, 0}, {cf.ng, cf.ngj, 2, 2});
  policy_f4 yPol = policy_f4({0, 0, 0, 0}, {cf.ngi, cf.ng, 2, 2});

  Kokkos::parallel_for(
      xPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
        lsR(i, j, k, l) = mt(mng + i, j, k, l);
        rsR(i, j, k, l) = mt(mnci + i, j, k, l);
      });
  Kokkos::parallel_for(
      yPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
        bsR(i, j, k, l) = mt(i, mng + j, k, l);
        tsR(i, j, k, l) = mt(i, mncj + j, k, l);
      });
  Kokkos::deep_copy(lsH, ls);
  Kokkos::deep_copy(rsH, rs);
  Kokkos::deep_copy(bsH, bs);
  Kokkos::deep_copy(tsH, ts);

  MPI_Request reqs[8];

  MPI_Isend(lsH.data(), cf.ng * cf.ngj * 4, MPI_DOUBLE, cf.xMinus, 0, cf.comm,
            &reqs[0]);
  MPI_Irecv(lrH.data(), cf.ng * cf.ngj * 4, MPI_DOUBLE, cf.xMinus, MPI_ANY_TAG,
            cf.comm, &reqs[1]);

  MPI_Isend(rsH.data(), cf.ng * cf.ngj * 4, MPI_DOUBLE, cf.xPlus, 0, cf.comm,
            &reqs[2]);
  MPI_Irecv(rrH.data(), cf.ng * cf.ngj * 4, MPI_DOUBLE, cf.xPlus, MPI_ANY_TAG,
            cf.comm, &reqs[3]);

  MPI_Isend(bsH.data(), cf.ngi * cf.ng * 4, MPI_DOUBLE, cf.yMinus, 0, cf.comm,
            &reqs[4]);
  MPI_Irecv(brH.data(), cf.ngi * cf.ng * 4, MPI_DOUBLE, cf.yMinus, MPI_ANY_TAG,
            cf.comm, &reqs[5]);

  MPI_Isend(tsH.data(), cf.ngi * cf.ng * 4, MPI_DOUBLE, cf.yPlus, 0, cf.comm,
            &reqs[6]);
  MPI_Irecv(trH.data(), cf.ngi * cf.ng * 4, MPI_DOUBLE, cf.yPlus, MPI_ANY_TAG,
            cf.comm, &reqs[7]);

  MPI_Waitall(8, reqs, MPI_STATUS_IGNORE);

  Kokkos::deep_copy(lr, lrH);
  Kokkos::deep_copy(rr, rrH);
  Kokkos::deep_copy(br, brH);
  Kokkos::deep_copy(tr, trH);

  Kokkos::parallel_for(
      xPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
        mt(i, j, k, l) = lrR(i, j, k, l);
        mt(mngi - mng + i, j, k, l) = rrR(i, j, k, l);
      });
  Kokkos::parallel_for(
      yPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
        mt(i, j, k, l) = brR(i, j, k, l);
        mt(i, mngj - mng + j, k, l) = trR(i, j, k, l);
      });
#endif

  if (cf.yMinus < 0) {
    Kokkos::parallel_for(
        "metricbcn", policy_f({0, 0}, {mng, mngi}),
        KOKKOS_LAMBDA(const int g, const int i) {
          mt(i, mng - g - 1, 0, 0) = mt(i, mng + g, 0, 0);
          mt(i, mng - g - 1, 0, 1) = mt(i, mng + g, 0, 1);
          mt(i, mng - g - 1, 1, 0) = mt(i, mng + g, 1, 0);
          mt(i, mng - g - 1, 1, 1) = mt(i, mng + g, 1, 1);
        });
  }
  if (cf.yPlus < 0) {
    Kokkos::parallel_for(
        "metricbct", policy_f({0, 0}, {mng, mngi}),
        KOKKOS_LAMBDA(const int g, const int i) {
          mt(i, mngj - 1 - g, 0, 0) = mt(i, mncj + g, 0, 0);
          mt(i, mngj - 1 - g, 0, 1) = mt(i, mncj + g, 0, 1);
          mt(i, mngj - 1 - g, 1, 0) = mt(i, mncj + g, 1, 0);
          mt(i, mngj - 1 - g, 1, 1) = mt(i, mncj + g, 1, 1);
        });
  }
  if (cf.xMinus < 0) {
    Kokkos::parallel_for(
        "metricbcl", policy_f({0, 0}, {mng, mngj}),
        KOKKOS_LAMBDA(const int g, const int j) {
          mt(mng - g - 1, j, 0, 0) = mt(mng + g, j, 0, 0);
          mt(mng - g - 1, j, 0, 1) = mt(mng + g, j, 0, 1);
          mt(mng - g - 1, j, 1, 0) = mt(mng + g, j, 1, 0);
          mt(mng - g - 1, j, 1, 1) = mt(mng + g, j, 1, 1);
        });
  }
  if (cf.xPlus < 0) {
    Kokkos::parallel_for(
        "metricbcr", policy_f({0, 0}, {mng, mngj}),
        KOKKOS_LAMBDA(const int g, const int j) {
          mt(mngi - 1 - g, j, 0, 0) = mt(mnci + g, j, 0, 0);
          mt(mngi - 1 - g, j, 0, 1) = mt(mnci + g, j, 0, 1);
          mt(mngi - 1 - g, j, 1, 0) = mt(mnci + g, j, 1, 0);
          mt(mngi - 1 - g, j, 1, 1) = mt(mnci + g, j, 1, 1);
        });
  }

  Kokkos::fence();
  timers["calcMetrics"].accumulate();
}

void gen2d_func::postSim() {}
