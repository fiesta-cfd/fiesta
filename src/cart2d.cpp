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

#include "cart2d.hpp"
#include "vtk.hpp"
#include "mpi.hpp"
#include <cassert>

std::map<string,int> varxIds;

cart2d_func::cart2d_func(struct inputConfig &cf_) : rk_func(cf_) {

  grid  = FS4D("coords", cf.ni, cf.nj, cf.nk, 3);          // Grid Coords
  var   = FS4D("var", cf.ngi, cf.ngj, cf.ngk, cf.nvt);     // Primary Variables
  tmp1  = FS4D("tmp1", cf.ngi, cf.ngj, cf.ngk, cf.nvt);    // Temporary Array
  dvar  = FS4D("dvar", cf.ngi, cf.ngj, cf.ngk, cf.nvt);    // RHS Output
  vel   = FS3D("vel", cf.ngi, cf.ngj, 2);                  // velocity
  p     = FS2D("p", cf.ngi, cf.ngj);                       // Pressure
  T     = FS2D("T", cf.ngi, cf.ngj);                       // Temperature
  rho   = FS2D("rho", cf.ngi, cf.ngj);                     // Total Density
  fluxx = FS2D("fluxx", cf.ngi, cf.ngj); // Advective Fluxes in X direction
  fluxy = FS2D("fluxy", cf.ngi, cf.ngj); // Advective Fluxes in Y direction

  varNames.push_back("X-Momentum");
  varNames.push_back("Y-Momentum");
  varNames.push_back("Energy");
  for (int v=0; v<cf.ns; ++v)
      varNames.push_back("Density " + cf.speciesName[v]);

  if (cf.visc == 1) {
    qx      = FS2D("qx", cf.ngi, cf.ngj); // Heat Fluxes X direction
    qy      = FS2D("qy", cf.ngi, cf.ngj); // Heat Fluxes Y direction
    stressx = FS4D("stressx", cf.ngi, cf.ngj, 2, 2); // stress on x faces
    stressy = FS4D("stressy", cf.ngi, cf.ngj, 2, 2); // stress on y faces
  }

  if (cf.ceq == 1) {
    gradRho = FS3D("gradRho", cf.ngi, cf.ngj, 4);    // Density Gradients
    varNames.push_back("C");
    varNames.push_back("C_hat");
    varNames.push_back("Tau_1");
    varNames.push_back("Tau_2");
    varNames.push_back("Unused");
  }

  if (cf.noise == 1) {
    noise = FS2D_I("noise", cf.ngi, cf.ngj);         // Noise Indicator
    varxNames.push_back("Noise");
  }

  assert(varNames.size() == cf.nvt);

  varxNames.push_back("X-Velocity");
  varxNames.push_back("Y-Velocity");
  varxNames.push_back("Pressure");
  varxNames.push_back("Temperature");
  varxNames.push_back("Density");

  varx  = FS4D("varx", cf.ngi, cf.ngj, cf.ngk, varxNames.size()); // Extra Vars

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

  // Create Simulation )timers
  timers["flux"] = fiestaTimer("Flux Calculation");
  timers["pressgrad"] = fiestaTimer("Pressure Gradient Calculation");
  timers["calcSecond"] = fiestaTimer("Secondary Variable Calculation");
  timers["solWrite"] = fiestaTimer("Solution Write Time");
  timers["resWrite"] = fiestaTimer("Restart Write Time");
  timers["statCheck"] = fiestaTimer("Status Check");
  timers["rk"] = fiestaTimer("Runge Stage Update");
  timers["halo"] = fiestaTimer("Halo Exchanges");
  timers["bc"] = fiestaTimer("Boundary Conditions");
  if (cf.gravity == 1) {
    timers["gravity"] = fiestaTimer("Gravity Term");
  }
  if (cf.visc == 1) {
    timers["stress"] = fiestaTimer("Stress Tensor Computation");
    timers["qflux"] = fiestaTimer("Heat Flux Calculation");
    timers["visc"] = fiestaTimer("Viscous Term Calculation");
  }
  if (cf.ceq == 1) {
    timers["ceq"] = fiestaTimer("C-Equation");
  }
  if (cf.noise == 1) {
    timers["noise"] = fiestaTimer("Noise Removal");
  }
};

void cart2d_func::preStep() {}

void cart2d_func::postStep() {

      policy_f ghost_pol = policy_f({0, 0}, {cf.ngi, cf.ngj});
  if (( (cf.write_freq >0) && (cf.t % cf.write_freq == 0) )||
      ( (cf.stat_freq  >0) && (cf.t % cf.stat_freq  == 0) )){

      // Calcualte Total Density and Pressure Fields
      timers["calcSecond"].reset();
      Kokkos::parallel_for(ghost_pol, calculateRhoPT2D(var, p, rho, T, cd));
      Kokkos::parallel_for(ghost_pol, computeVelocity2D(var, rho, vel));
      Kokkos::parallel_for(ghost_pol, copyExtraVars2D(varx, vel, p, rho, T));
      Kokkos::fence();
      timers["calcSecond"].accumulate();
    }

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

    // double etat = 5.0e-3;
    // double dh = 5.0e-4;
    // double doff = 0.1;

    policy_f noise_pol = policy_f({0, 0}, {M, N});
    policy_f cell_pol =
        policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});

    if (cf.ceq == 1) {
      double maxCh, maxCh_recv;
      Kokkos::parallel_reduce(cell_pol, maxCvar2D(var, 1, cd),
                              Kokkos::Max<double>(maxCh));
#ifndef NOMPI
      MPI_Allreduce(&maxCh, &maxCh_recv, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
      maxCh = maxCh_recv;
#endif
      coff = cf.n_coff * maxCh;

    } else {
      coff = 0.0;
    }

    timers["noise"].reset();
    for (int v = 0; v < 2; ++v) {
      // int v = 2;
      Kokkos::parallel_for(noise_pol,
                           detectNoise2D(var, noise, cf.n_dh, coff, cd, v));
      for (int tau = 0; tau < cf.n_nt; ++tau) {
        Kokkos::parallel_for(cell_pol,
                             removeNoise2D(dvar, var, noise, cf.n_eta, cd, v));
        Kokkos::parallel_for(cell_pol, updateNoise2D(dvar, var, v));
      }
    }
    Kokkos::fence();
    timers["noise"].accumulate();
  } // end noise
}

void cart2d_func::preSim() {
}
void cart2d_func::postSim() {}

void cart2d_func::compute() {

  policy_f ghost_pol = policy_f({0, 0}, {cf.ngi, cf.ngj});
  policy_f cell_pol =
      policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});
  policy_f face_pol =
      policy_f({cf.ng - 1, cf.ng - 1}, {cf.ngi - cf.ng, cf.ngj - cf.ng});

  // Calcualte Total Density and Pressure Fields
  timers["calcSecond"].reset();
  Kokkos::parallel_for(ghost_pol, calculateRhoPT2D(var, p, rho, T, cd));
  Kokkos::parallel_for(ghost_pol, computeVelocity2D(var, rho, vel));
  Kokkos::fence();
  timers["calcSecond"].accumulate();

  // Calculate and apply weno fluxes for each variable
  for (int v = 0; v < cf.nv; ++v) {
    timers["flux"].reset();
    if (cf.scheme == 3) {
      Kokkos::parallel_for(
          face_pol, computeFluxQuick2D(var, p, rho, fluxx, fluxy, cd, v));
    } else if (cf.scheme == 2) {
      Kokkos::parallel_for(
          face_pol, computeFluxCentered2D(var, p, rho, fluxx, fluxy, cd, v));
    } else {
      Kokkos::parallel_for(
          face_pol, computeFluxWeno2D(var, p, rho, vel, fluxx, fluxy, cd, v));
    }
    Kokkos::fence();
    Kokkos::parallel_for(cell_pol, advect2D(dvar, fluxx, fluxy, cd, v));
    Kokkos::fence();
    timers["flux"].accumulate();
  }

  // Apply Pressure Gradient Term
  timers["pressgrad"].reset();
  Kokkos::parallel_for(cell_pol, applyPressureGradient2D(dvar, p, cd));
  Kokkos::fence();
  timers["pressgrad"].accumulate();

  if (cf.gravity == 1) {
    timers["gravity"].reset();
    Kokkos::parallel_for(cell_pol, computeBuoyancy(dvar, var, rho));
    Kokkos::fence();
    timers["gravity"].accumulate();
  }

  if (cf.visc == 1) {
    timers["stress"].reset();
    Kokkos::parallel_for(
        face_pol, calculateStressTensor2dv(var, rho, vel, stressx, stressy, cd));
    Kokkos::fence();
    timers["stress"].accumulate();

    timers["qflux"].reset();
    Kokkos::parallel_for(face_pol,
                         calculateHeatFlux2dv(var, rho, T, qx, qy, cd));
    Kokkos::fence();
    timers["qflux"].accumulate();

    timers["visc"].reset();
    Kokkos::parallel_for(cell_pol, applyViscousTerm2dv(dvar, var, rho, vel, stressx,
                                                       stressy, qx, qy, cd));
    Kokkos::fence();
    timers["visc"].accumulate();
  }

  if (cf.ceq == 1) {
    double mu;
    double maxS,    maxS_recv;
    double maxG,    maxG_recv;
    double maxGh,   maxGh_recv;
    double maxTau1, maxTau1_recv;
    double maxTau2, maxTau2_recv;
    double maxC,    maxC_recv;
    double maxCh,   maxCh_recv;
    double maxT1,   maxT1_recv;
    double maxT2,   maxT2_recv;

    timers["ceq"].reset();
    Kokkos::parallel_reduce(cell_pol, maxWaveSpeed2D(var, p, rho, cd),
                            Kokkos::Max<double>(maxS));
    Kokkos::parallel_reduce(cell_pol, maxGradRho2D(gradRho, 1),
                            Kokkos::Max<double>(maxG));
    Kokkos::parallel_reduce(cell_pol, maxGradRho2D(gradRho, 1),
                            Kokkos::Max<double>(maxGh));
    Kokkos::parallel_reduce(cell_pol, maxGradRho2D(gradRho, 1),
                            Kokkos::Max<double>(maxTau1));
    Kokkos::parallel_reduce(cell_pol, maxGradRho2D(gradRho, 1),
                            Kokkos::Max<double>(maxTau2));
    Kokkos::parallel_reduce(cell_pol, maxCvar2D(var, 0, cd),
                            Kokkos::Max<double>(maxC));
    Kokkos::parallel_reduce(cell_pol, maxCvar2D(var, 1, cd),
                            Kokkos::Max<double>(maxCh));
    Kokkos::parallel_reduce(cell_pol, maxCvar2D(var, 2, cd),
                            Kokkos::Max<double>(maxT1));
    Kokkos::parallel_reduce(cell_pol, maxCvar2D(var, 3, cd),
                            Kokkos::Max<double>(maxT2));

#ifndef NOMPI
    MPI_Allreduce(&maxS, &maxS_recv, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxG, &maxG_recv, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxGh, &maxGh_recv, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxTau1, &maxTau1_recv, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxTau2, &maxTau2_recv, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxC, &maxC_recv, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxCh, &maxCh_recv, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxT1, &maxT1_recv, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxT2, &maxT2_recv, 1, MPI_DOUBLE, MPI_MAX, cf.comm);

    maxS    = maxS_recv;
    maxG    = maxG_recv;
    maxGh   = maxGh_recv;
    maxTau1 = maxTau1_recv;
    maxTau2 = maxTau2_recv;
    maxC    = maxC_recv;
    maxCh   = maxCh_recv;
    maxT1   = maxT1_recv;
    maxT2   = maxT2_recv;
#endif

    if (maxG < 1e-6)
      maxG = 1e-6;
    if (maxGh < 1e-6)
      maxGh = 1e-6;
    if (maxTau1 < 1e-6)
      maxTau1 = 1e-6;
    if (maxTau2 < 1e-6)
      maxTau2 = 1e-6;
    if (maxC < 1e-6)
      maxC = 1e-6;
    if (maxCh < 1e-6)
      maxCh = 1e-6;
    if (maxT1 < 1e-6)
      maxT1 = 1e-6;
    if (maxT2 < 1e-6)
      maxT2 = 1e-6;

    mu = maxT1;
    if (maxT2 > mu)
      mu = maxT2;

    // printf("### DEBUG ### %f : %f : %f : %f :
    // %f\n",maxC,maxCh,maxTau1,maxTau2,mu);

    Kokkos::parallel_for(cell_pol, calculateRhoGrad2D(var, rho, gradRho, cd));

    Kokkos::parallel_for(cell_pol,
                         updateCeq2D(dvar, var, gradRho, maxS, cd, cf.kap,
                                     cf.eps, maxG, maxGh, maxTau1, maxTau2));

    if (cf.t >= cf.st){
      Kokkos::parallel_for(cell_pol, applyCeq2D(dvar, var, rho, cf.beta, cf.betae,
                                                cf.alpha, maxC, maxCh, mu, cd));
    }

    Kokkos::fence();
    timers["ceq"].accumulate();
  }
}
