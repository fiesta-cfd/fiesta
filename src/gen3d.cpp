#include "fiesta.hpp"
#include "input.hpp"
#ifndef NOMPI
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

gen3d_func::gen3d_func(struct inputConfig &cf_, Kokkos::View<double *> &cd_)
    : rk_func(cf_, cd_) {

  // Allocate all device arrays here
  grid    = FS4D("coords",  cf.ni,  cf.nj,  cf.nk,  3);
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
  cd = mcd;

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
  timers["flux"]        = fiestaTimer("Flux Calculation");
  timers["pressgrad"]   = fiestaTimer("Pressure Gradient Calculation");
  timers["calcSecond"]  = fiestaTimer("Secondary Variable Calculation");
  timers["solWrite"]    = fiestaTimer("Solution Write Time");
  timers["resWrite"]    = fiestaTimer("Restart Write Time");
  timers["statCheck"]   = fiestaTimer("Status Check");
  timers["calcMetrics"] = fiestaTimer("Metric Computations");
};

void gen3d_func::preStep() {}

void gen3d_func::postStep() {
  
  // MDRange Policy for all cells including ghost cells
  policy_f3 ghost_pol = policy_f3({0, 0, 0}, {cf.ngi, cf.ngj, cf.ngk});

  // Copy secondary variables to extra variables array
  if (cf.write_freq >0)
    if (cf.t % cf.write_freq == 0){

      timers["calcSecond"].reset();
      Kokkos::parallel_for(ghost_pol, calculateRhoPT3D(var, p, rho, T, cd));
      Kokkos::parallel_for(ghost_pol, computeVelocity3D(var, rho, tvel));
      Kokkos::fence();
      timers["calcSecond"].accumulate();

      Kokkos::parallel_for( "CopyVarx", ghost_pol,
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          varx(i,j,k,0) = tvel(i,j,k,0);
          varx(i,j,k,1) = tvel(i,j,k,1);
          varx(i,j,k,2) = tvel(i,j,k,2);
          varx(i,j,k,3) =    p(i,j,k);
          varx(i,j,k,4) =    T(i,j,k);
          varx(i,j,k,5) =  rho(i,j,k);
      });
    }
}
void gen3d_func::preSim() {
  policy_f3 cell_pol3 = policy_f3(
      {cf.ng, cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng, cf.ngk - cf.ng});

  // compute metrics
  timers["calcMetrics"].reset();
  Kokkos::parallel_for(cell_pol3, computeMetrics3D(metrics, grid));

  FS5D &mt = metrics;
  inputConfig &c = cf;

#ifndef NOMPI
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
        lsR(i, j, k, m, n) = mt(c.ng + i, j, k, m, n);
        rsR(i, j, k, m, n) = mt(c.nci + i, j, k, m, n);
      });
  Kokkos::parallel_for(
      yPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int m,
                          const int n) {
        bsR(i, j, k, m, n) = mt(i, c.ng + j, k, m, n);
        tsR(i, j, k, m, n) = mt(i, c.ncj + j, k, m, n);
      });
  Kokkos::parallel_for(
      zPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int m,
                          const int n) {
        hsR(i, j, k, m, n) = mt(i, j, c.ng + k, m, n);
        fsR(i, j, k, m, n) = mt(i, j, c.nck + k, m, n);
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
        mt(c.ngi - c.ng + i, j, k, m, n) = rrR(i, j, k, m, n);
      });
  Kokkos::parallel_for(
      yPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int m,
                          const int n) {
        mt(i, j, k, m, n) = brR(i, j, k, m, n);
        mt(i, c.ngj - c.ng + j, k, m, n) = trR(i, j, k, m, n);
      });
  Kokkos::parallel_for(
      zPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int m,
                          const int n) {
        mt(i, j, k, m, n) = hrR(i, j, k, m, n);
        mt(i, j, c.ngk - c.ng + k, m, n) = frR(i, j, k, m, n);
      });

#endif

  if (cf.xMinus < 0) {
    Kokkos::parallel_for(
        "metricbcl", policy_f3({0, 0, 0}, {cf.ng, cf.ngj, cf.ngk}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
              mt(c.ng - i - 1, j, k, m, n) = mt(c.ng + i, j, k, m, n);
        });
  }
  if (cf.xPlus < 0) {
    Kokkos::parallel_for(
        "metricbcr", policy_f3({0, 0, 0}, {cf.ng, cf.ngj, cf.ngk}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
              mt(c.ngi - 1 - i, j, k, m, n) = mt(c.nci + i, j, k, m, n);
        });
  }
  if (cf.yMinus < 0) {
    Kokkos::parallel_for(
        "metricbcb", policy_f3({0, 0, 0}, {cf.ngi, cf.ng, cf.ngk}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
              mt(i, c.ng - j - 1, k, m, n) = mt(i, c.ng + j, k, m, n);
        });
  }
  if (cf.yPlus < 0) {
    Kokkos::parallel_for(
        "metricbct", policy_f3({0, 0, 0}, {cf.ngi, cf.ng, cf.ngk}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
              mt(i, c.ngj - 1 - j, k, m, n) = mt(i, c.ncj + j, k, m, n);
        });
  }
  if (cf.zMinus < 0) {
    Kokkos::parallel_for(
        "metricbch", policy_f3({0, 0, 0}, {cf.ngi, cf.ngj, cf.ng}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
              mt(i, j, c.ng - k - 1, m, n) = mt(i, j, c.ng + k, m, n);
        });
  }
  if (cf.zPlus < 0) {
    Kokkos::parallel_for(
        "metricbcf", policy_f3({0, 0, 0}, {cf.ngi, cf.ngj, cf.ng}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
              mt(i, j, c.ngk - 1 - k, m, n) = mt(i, j, c.nck + k, m, n);
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
        calculateWenoFluxesG(var, p, rho, tvel, fluxx, fluxy, fluxz, cd, v));
    Kokkos::parallel_for(cell_pol, advect3D(dvar, fluxx, fluxy, fluxz, v));
  }
  Kokkos::fence();
  timers["flux"].accumulate();

  timers["pressgrad"].reset();
  Kokkos::parallel_for(cell_pol,
                       applyGenPressureGradient3D(dvar, metrics, p, cd));
  Kokkos::fence();
  timers["pressgrad"].accumulate();
}
