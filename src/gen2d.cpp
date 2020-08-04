#include "gen2d.hpp"

gen2d_func::gen2d_func(struct inputConfig &cf_, Kokkos::View<double *> &cd_)
    : rk_func(cf_, cd_) {
  /* Generalized 2D module constructor */

  // Data Structures (see fiesta.hpp to decode typenames)
  grid = FS4D("grid", cf.ni, cf.nj, cf.nk, 3);       // grid coordinates
  metrics = FS4D("metrics", cf.ngi, cf.ngj, 2, 2);   // Jacobian metrics
  var = FS4D("var", cf.ngi, cf.ngj, cf.ngk, cf.nvt); // Primary Variable Array
  tmp1 = FS4D("tmp1", cf.ngi, cf.ngj, cf.ngk,
              cf.nvt); // Temporary Variable (for Runge-Kutta)
  dvar = FS4D("dvar", cf.ngi, cf.ngj, cf.ngk, cf.nvt); // RHS Output
  tvel = FS3D("vel", cf.ngi, cf.ngj, 2); // Transformed velocity components
  p = FS2D("p", cf.ngi, cf.ngj);         // Pressure
  T = FS2D("T", cf.ngi, cf.ngj);         // Temperature
  rho = FS2D("rho", cf.ngi, cf.ngj);     // Total Density
  fluxx = FS2D("fluxx", cf.ngi, cf.ngj); // Advective Fluxes in X direction
  fluxy = FS2D("fluxy", cf.ngi, cf.ngj); // Advective Fluxes in Y direction
  if (cf.noise == 1) {
    noise = FS2D_I("noise", cf.ngi, cf.ngj); // Noise indicator array
  }
#ifdef NOMPI
  if (cf.particle == 1) {
    particles = FSP2D("particles", cf.p_np); // 2D particle view
    particlesH = Kokkos::create_mirror_view(particles);
  }
#endif
  cd = mcd;

  // Set up timer dictionaries
  timers["flux"] = fiestaTimer("Flux Calculation");
  timers["pressgrad"] = fiestaTimer("Pressure Gradient Calculation");
  timers["calcSecond"] = fiestaTimer("Secondary Variable Calculation");
  timers["solWrite"] = fiestaTimer("Solution Write Time");
  timers["resWrite"] = fiestaTimer("Restart Write Time");
  timers["statCheck"] = fiestaTimer("Status Check");
  timers["rk"] = fiestaTimer("Runge Stage Update");
  timers["calcMatrics"] = fiestaTimer("Metric Calculations");
  if (cf.noise == 1) {
    timers["noise"] = fiestaTimer("Noise Removal");
  }
  if (cf.particle == 1) {
    timers["padvect"] = fiestaTimer("Particle Advection");
    timers["pwrite"] = fiestaTimer("Particle Write");
    timers["psetup"] = fiestaTimer("Particle Setup");
  }

  ghostPol = policy_f({0, 0}, {cf.ngi, cf.ngj});
  cellPol = policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});
  facePol = policy_f({cf.ng - 1, cf.ng - 1}, {cf.ngi - cf.ng, cf.ngj - cf.ng});
};

void gen2d_func::compute() {
  /* main compute function for generalized 2d module */

  // Calcualte Total Density and Pressure Fields
  timers["calcSecond"].reset();
  Kokkos::parallel_for(ghostPol, calculateRhoPT2D(var, p, rho, T, cd));
  Kokkos::parallel_for(ghostPol, computeGenVelocity2D(var, metrics, rho, tvel));
  Kokkos::fence();
  timers["calcSecond"].accumulate();

  // Calculate and apply weno fluxes for each variable
  timers["flux"].reset();
  for (int v = 0; v < cf.nv; ++v) {
    if (cf.scheme == 3) {
      Kokkos::parallel_for(
          facePol, computeFluxQuick2D(var, p, rho, fluxx, fluxy, cd, v));
    } else if (cf.scheme == 2) {
      Kokkos::parallel_for(
          facePol, computeFluxCentered2D(var, p, rho, fluxx, fluxy, cd, v));
    } else {
      Kokkos::parallel_for(
          facePol, computeFluxWeno2D(var, p, rho, tvel, fluxx, fluxy, cd, v));
    }
    Kokkos::parallel_for(cellPol, advect2D(dvar, fluxx, fluxy, cd, v));
  }
  Kokkos::fence();
  timers["flux"].accumulate();

  // Apply Pressure Gradient Term
  timers["pressgrad"].reset();
  Kokkos::parallel_for(cellPol,
                       applyGenPressureGradient2D(dvar, metrics, p, cd));
  Kokkos::fence();
  timers["pressgrad"].accumulate();
}

void gen2d_func::preStep() {
  /* Pre time step hook executed before the first runge kutta stage */
}

void gen2d_func::postStep() {
  /* Post time step hook executed after the last runge kutta stage */

  if (cf.noise == 1) {
    int M = 0;
    int N = 0;
    double coff;

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
                           detectNoise2D(var, noise, cf.n_dh, coff, cd, v));
      for (int tau = 0; tau < cf.n_nt; ++tau) {
        Kokkos::parallel_for(cellPol,
                             removeNoise2D(dvar, var, noise, cf.n_eta, cd, v));
        Kokkos::parallel_for(cellPol, updateNoise2D(dvar, var, v));
      }
    }
    Kokkos::fence();
    timers["noise"].accumulate();
  } // end noise

#ifdef NOMPI
  if (cf.particle == 1) {
    // write particle data
    if (cf.write_freq > 0) {
      if ((cf.t) % cf.write_freq == 0) {
        timers["pwrite"].reset();
        Kokkos::deep_copy(particlesH, particles);
        // writeParticles(cf,particlesH);
        Kokkos::fence();
        timers["pwrite"].accumulate();
      }
    } // end particle write

    // execution policy for all cells including ghost cells
    policy_f ghost_pol = policy_f({0, 0}, {cf.ngi, cf.ngj});

    // Calcualte Total Density and Pressure Fields
    timers["calcSecond"].reset();
    Kokkos::parallel_for(ghost_pol, calculateRhoPT2D(var, p, rho, T, cd));
    Kokkos::fence();
    timers["calcSecond"].accumulate();

    // advect particles
    timers["padvect"].reset();
    Kokkos::parallel_for(
        cf.p_np, advectParticles2D(var, rho, grid, particles, cf.dt, cf.ng));
    Kokkos::fence();
    timers["padvect"].accumulate();
  }
#endif
}
void gen2d_func::preSim() {

  // execution policy for all cells including ghost cells
  policy_f cell_pol =
      policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});

  // compute metrics
  timers["calcMetrics"].reset();
  Kokkos::parallel_for(cell_pol, computeMetrics2D(metrics, grid));

  FS4D &mt = metrics;
  inputConfig &c = cf;
#ifndef NOMPI
  // mpi exchange of metrics
  ls = Kokkos::View<double ****, FS_LAYOUT>("leftSend", cf.ng, cf.ngj, 2, 2);
  lr = Kokkos::View<double ****, FS_LAYOUT>("leftRecv", cf.ng, cf.ngj, 2, 2);
  rs = Kokkos::View<double ****, FS_LAYOUT>("rightSend", cf.ng, cf.ngj, 2, 2);
  rr = Kokkos::View<double ****, FS_LAYOUT>("rightRecv", cf.ng, cf.ngj, 2, 2);
  bs = Kokkos::View<double ****, FS_LAYOUT>("bottomSend", cf.ngi, cf.ng, 2, 2);
  br = Kokkos::View<double ****, FS_LAYOUT>("bottomRecv", cf.ngi, cf.ng, 2, 2);
  ts = Kokkos::View<double ****, FS_LAYOUT>("topSend", cf.ngi, cf.ng, 2, 2);
  tr = Kokkos::View<double ****, FS_LAYOUT>("topRecv", cf.ngi, cf.ng, 2, 2);

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
        lsR(i, j, k, l) = mt(cf.ng + i, j, k, l);
        rsR(i, j, k, l) = mt(cf.nci + i, j, k, l);
      });
  Kokkos::parallel_for(
      yPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
        bsR(i, j, k, l) = mt(i, cf.ng + j, k, l);
        tsR(i, j, k, l) = mt(i, cf.ncj + j, k, l);
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
        mt(cf.ngi - cf.ng + i, j, k, l) = rrR(i, j, k, l);
      });
  Kokkos::parallel_for(
      yPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
        mt(i, j, k, l) = brR(i, j, k, l);
        mt(i, cf.ngj - cf.ng + j, k, l) = trR(i, j, k, l);
      });
#endif

  if (cf.yMinus < 0) {
    Kokkos::parallel_for(
        "metricbcn", policy_f({0, 0}, {cf.ng, cf.ngi}),
        KOKKOS_LAMBDA(const int g, const int i) {
          mt(i, c.ng - g - 1, 0, 0) = mt(i, c.ng + g, 0, 0);
          mt(i, c.ng - g - 1, 0, 1) = mt(i, c.ng + g, 0, 1);
          mt(i, c.ng - g - 1, 1, 0) = mt(i, c.ng + g, 1, 0);
          mt(i, c.ng - g - 1, 1, 1) = mt(i, c.ng + g, 1, 1);
        });
  }
  if (cf.yPlus < 0) {
    Kokkos::parallel_for(
        "metricbct", policy_f({0, 0}, {cf.ng, cf.ngi}),
        KOKKOS_LAMBDA(const int g, const int i) {
          mt(i, c.ngj - 1 - g, 0, 0) = mt(i, c.ncj + g, 0, 0);
          mt(i, c.ngj - 1 - g, 0, 1) = mt(i, c.ncj + g, 0, 1);
          mt(i, c.ngj - 1 - g, 1, 0) = mt(i, c.ncj + g, 1, 0);
          mt(i, c.ngj - 1 - g, 1, 1) = mt(i, c.ncj + g, 1, 1);
        });
  }
  if (cf.xMinus < 0) {
    Kokkos::parallel_for(
        "metricbcl", policy_f({0, 0}, {cf.ng, cf.ngj}),
        KOKKOS_LAMBDA(const int g, const int j) {
          mt(c.ng - g - 1, j, 0, 0) = mt(c.ng + g, j, 0, 0);
          mt(c.ng - g - 1, j, 0, 1) = mt(c.ng + g, j, 0, 1);
          mt(c.ng - g - 1, j, 1, 0) = mt(c.ng + g, j, 1, 0);
          mt(c.ng - g - 1, j, 1, 1) = mt(c.ng + g, j, 1, 1);
        });
  }
  if (cf.xPlus < 0) {
    Kokkos::parallel_for(
        "metricbcr", policy_f({0, 0}, {cf.ng, cf.ngj}),
        KOKKOS_LAMBDA(const int g, const int j) {
          mt(c.ngi - 1 - g, j, 0, 0) = mt(c.nci + g, j, 0, 0);
          mt(c.ngi - 1 - g, j, 0, 1) = mt(c.nci + g, j, 0, 1);
          mt(c.ngi - 1 - g, j, 1, 0) = mt(c.nci + g, j, 1, 0);
          mt(c.ngi - 1 - g, j, 1, 1) = mt(c.nci + g, j, 1, 1);
        });
  }

  Kokkos::fence();
  timers["calcMetrics"].accumulate();

  if (cf.particle == 1) {
    timers["psetup"].reset();
    Kokkos::View<particleStruct2D *>::HostMirror particlesH =
        Kokkos::create_mirror_view(particles);

    double ymax = 6.0;
    double xmax = 2.0;
    double dpx = xmax / ((double)(cf.p_np - 1));
    double dpy = 9.0 / ((double)(cf.p_np - 1));
    for (int p = 0; p < cf.p_np; ++p) {
      particlesH(p).state = 1;
      particlesH(p).x = 1.0;
      particlesH(p).y = 0.5 + p * dpy;
    }

    // find initial cell id
    policy_f grid_pol = policy_f({0, 0}, {cf.nci, cf.ncj});
    for (int p = 0; p < cf.p_np; ++p) {
      Kokkos::parallel_for(grid_pol,
                           findInitialCell2D(grid, particles, p, cf.ng));
    }
    Kokkos::fence();
    timers["psetup"].accumulate();

    Kokkos::deep_copy(particles, particlesH);

    if (cf.write_freq > 0) {
      timers["pwrite"].reset();
      Kokkos::deep_copy(particlesH, particles);
      // writeParticles(cf,particlesH);
      Kokkos::fence();
      timers["pwrite"].accumulate();
    } // end initial write
  }
}
void gen2d_func::postSim() {}
