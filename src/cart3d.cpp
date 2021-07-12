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
#ifndef NOMPI
#include "mpi.hpp"
#include "mpi.h"
#endif
#include "Kokkos_Core.hpp"
#include "bc.hpp"
#include "cart3d.hpp"
#include "debug.hpp"
#include "secondary.hpp"
#include "velocity.hpp"
#include "flux.hpp"
#include "advect.hpp"
#include "presgrad.hpp"
#include "block.hpp"
#include "ceq.hpp"
#include "noise.hpp"
#include "log2.hpp"
#include "buoyancy.hpp"

cart3d_func::cart3d_func(struct inputConfig &cf_) : rk_func(cf_) {
  // Allocate all device variables here

  size_t memEstimate = 3*cf.nvt+6+3+6;
  if (cf.visc==1) memEstimate += 3+27;
  if (cf.ceq==1) memEstimate += 5+3+3;

  memEstimate *= (cf.ngi*cf.ngj*cf.ngk);
  memEstimate += 3*(cf.ni+cf.nj+cf.nk);
  memEstimate *= 8;  // 8 bytes per double
  double memEstimateMB = memEstimate/(1048576.0);

  cf.log->message("Minimum device memory estimate: {:.2f}",memEstimateMB);
  cf.log->message("Minimum device memory estimate: {}",memEstimate);

  grid    = FS4D("coords",    cf.ni,  cf.nj,  cf.nk,  3);      // Grid Coords

  var     = FS4D("var",       cf.ngi, cf.ngj, cf.ngk, cf.nvt); // Primary Vars
  tmp1    = FS4D( "tmp1",     cf.ngi, cf.ngj, cf.ngk, cf.nvt); // Temp Vars
  dvar    = FS4D("dvar",      cf.ngi, cf.ngj, cf.ngk, cf.nvt); // RHS Output

  p       = FS3D("p",         cf.ngi, cf.ngj, cf.ngk);         // Pressure
  T       = FS3D("T",         cf.ngi, cf.ngj, cf.ngk);         // Temperature
  rho     = FS3D("rho",       cf.ngi, cf.ngj, cf.ngk);         // Total Density
  vel     = FS4D("vel",       cf.ngi, cf.ngj, cf.ngk,3);         // Total Density

  fluxx   = FS3D( "fluxx",    cf.ngi, cf.ngj, cf.ngk);    // Advective Fluxes X
  fluxy   = FS3D( "fluxy",    cf.ngi, cf.ngj, cf.ngk);    // Advective Fluxes Y
  fluxz   = FS3D( "fluxz",    cf.ngi, cf.ngj, cf.ngk);    // Advective Fluxes Z

  if (cf.visc == 1) {
    qx      = FS3D( "qx",       cf.ngi, cf.ngj, cf.ngk);         // Heat Flux X
    qy      = FS3D( "qy",       cf.ngi, cf.ngj, cf.ngk);         // Heat Flux Y
    qz      = FS3D( "qz",       cf.ngi, cf.ngj, cf.ngk);         // Heat Flux Z
    stressx = FS5D( "stressx",  cf.ngi, cf.ngj, cf.ngk, 3, 3); // stress tensor X
    stressy = FS5D( "stressy",  cf.ngi, cf.ngj, cf.ngk, 3, 3); // stress tensor Y
    stressz = FS5D( "stressz",  cf.ngi, cf.ngj, cf.ngk, 3, 3); // stress tensor Z
  }
  if (cf.ceq) {
    gradRho = FS4D( "gradRho",  cf.ngi, cf.ngj, cf.ngk, 5);    // Density Gradien
    cFlux   = FS4D("cFlux",     cf.ngi, cf.ngj, cf.ngk, 3);    // 
    mFlux   = FS6D("mFlux", 3,3,cf.ngi, cf.ngj, cf.ngk, 3);    //
  }
  if (cf.noise) {
    noise = FS3D_I("noise", cf.ngi, cf.ngj, cf.ngk);
  }

  // Primary Variable Names
  varNames.push_back("X-Momentum");
  varNames.push_back("Y-Momentum");
  varNames.push_back("Z-Momentum");
  varNames.push_back("Energy");
  for (int v=0; v<cf.ns; ++v)
    varNames.push_back(cf.speciesName[v] + " Density");
  if (cf.ceq == 1) {
    varNames.push_back("C");
    varNames.push_back("C_hat");
    varNames.push_back("C_1");
    varNames.push_back("C_2");
    varNames.push_back("C_3");
  }

  assert(varNames.size()==cf.nvt);

  // Secondary Variable Names
  varxNames.push_back("X-Velocity");
  varxNames.push_back("Y-Velocity");
  varxNames.push_back("Z-Velocity");
  varxNames.push_back("Pressure");
  varxNames.push_back("Temperature");
  varxNames.push_back("Total Density");
  if (cf.noise){
    varxNames.push_back("Noise_c");
    varxNames.push_back("Noise_I");
    varxNames.push_back("Noise_d");
  }
  if (cf.ceq)
    varxNames.push_back("dc");
  Fiesta::Log::debug("varxName.size()={}",varxNames.size());

  // Create Secondary Variable Array
  varx = FS4D("varx",cf.ngi,cf.ngj,cf.ngk,varxNames.size());


  // Create Timers
  timers["flux"] = fiestaTimer("Flux Calculation");
  timers["pressgrad"] = fiestaTimer("Pressure Gradient Calculation");
  timers["calcSecond"] = fiestaTimer("Secondary Variable Calculation");
  timers["solWrite"] = fiestaTimer("Solution Write Time");
  timers["resWrite"] = fiestaTimer("Restart Write Time");
  timers["statCheck"] = fiestaTimer("Status Check");
  timers["rk"] = fiestaTimer("Runge Stage Update");
  timers["halo"] = fiestaTimer("Halo Exchanges");
  timers["bc"] = fiestaTimer("Boundary Conditions");
  if (cf.visc == 1) {
    timers["stress"] = fiestaTimer("Stress Tensor Computation");
    timers["qflux"] = fiestaTimer("Heat Flux Calculation");
    timers["visc"] = fiestaTimer("Viscous Term Calculation");
  }
  if (cf.ceq) {
    timers["ceq"] = fiestaTimer("C-Equation");
  }
  if (cf.noise) {
    timers["noise"] = fiestaTimer("Noise Filter");
  }
  if (cf.buoyancy) {
    timers["buoyancy"] = fiestaTimer("Buoyancy Term");
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
};

void cart3d_func::preStep() {}

void cart3d_func::postStep() {
  bool varsxNeeded=false;

  //for (auto& block : cf.ioblocks)
  //  if( (block.frq() > 0) && (cf.t % block.frq() == 0) ) varsxNeeded=true;

  if ( (cf.write_freq >0) && (cf.t % cf.write_freq == 0) ) varsxNeeded=true;
  if ( (cf.stat_freq  >0) && (cf.t % cf.stat_freq  == 0) ) varsxNeeded=true;

  if (varsxNeeded){
    timers["calcSecond"].reset();
    policy_f3 ghost_pol = policy_f3({0, 0, 0}, {cf.ngi, cf.ngj, cf.ngk});
    Kokkos::parallel_for(ghost_pol, calculateRhoPT3D(var, p, rho, T, cd));
    Kokkos::parallel_for(ghost_pol, computeVelocity3D(var, rho, vel));
    Kokkos::parallel_for(ghost_pol, copyExtraVars3D(varx, vel, p, rho, T));
    Kokkos::fence();
    timers["calcSecond"].accumulate();
  }

  if (cf.noise) {
    timers["noise"].reset();
    int M = 0;
    int N = 0;
    int L = 0;
    double coff;

    if ((cf.nci - 1) % 2 == 0)
      M = (cf.nci - 1) / 2;
    else
      M = cf.nci / 2;

    if ((cf.ncj - 1) % 2 == 0)
      N = (cf.ncj - 1) / 2;
    else
      N = cf.ncj / 2;

    if ((cf.nck - 1) % 2 == 0)
      L = (cf.nck - 1) / 2;
    else
      L = cf.nck / 2;

    policy_f3 noise_pol = policy_f3({0, 0, 0}, {M, N, L});
    policy_f3 cell_pol = policy_f3({cf.ng, cf.ng, cf.ng}, {cf.ngi-cf.ng, cf.ngj-cf.ng, cf.ngk-cf.ng});

    if (cf.ceq == 1) {
      double maxCh=0;
      Kokkos::parallel_reduce(cell_pol, maxCvar3D(var, 6), Kokkos::Max<double>(maxCh));
      #ifndef NOMPI
        MPI_Allreduce(&maxCh, &maxCh, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
      #endif
      coff = cf.n_coff * maxCh;
    } else {
      coff = 0.0;
    }

    for (int v = 0; v < 3; ++v) {
      Kokkos::parallel_for(noise_pol, detectNoise3D(var, varx, noise, cf.n_dh, coff, cd, v));
      for (int tau = 0; tau < cf.n_nt; ++tau) {
        Kokkos::parallel_for(cell_pol, removeNoise3D(dvar, var, varx, noise, cf.n_eta, cd, v));
        Kokkos::parallel_for(cell_pol, updateNoise3D(dvar, var, v));
      }
    }
    Kokkos::fence();
    timers["noise"].accumulate();
  } // end noise
}

void cart3d_func::preSim() {
  policy_f3 ghost_pol = policy_f3({0, 0, 0}, {cf.ngi, cf.ngj, cf.ngk});

  timers["calcSecond"].reset();
  Kokkos::parallel_for(ghost_pol, calculateRhoPT3D(var, p, rho, T, cd));
  Kokkos::parallel_for(ghost_pol, computeVelocity3D(var, rho, vel));
  Kokkos::parallel_for(ghost_pol, copyExtraVars3D(varx, vel, p, rho, T));
  Kokkos::fence();
  timers["calcSecond"].accumulate();
}

void cart3d_func::postSim() {}

void cart3d_func::compute() {
  // create range policies
  policy_f3 ghost_pol = policy_f3({0, 0, 0}, {cf.ngi, cf.ngj, cf.ngk});
  policy_f3 cell_pol = policy_f3(
      {cf.ng, cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng, cf.ngk - cf.ng});
  policy_f3 weno_pol =
      policy_f3({cf.ng - 1, cf.ng - 1, cf.ng - 1},
                {cf.ngi - cf.ng, cf.ngj - cf.ng, cf.ngk - cf.ng});

  double maxS;

  double maxC;
  double maxCh;
  double maxC1;
  double maxC2;
  double maxC3;
  double mu, alpha, beta, betae;

  /**** WENO ****/
  timers["calcSecond"].reset();
  Kokkos::parallel_for(ghost_pol, calculateRhoPT3D(var, p, rho, T, cd));
  Kokkos::parallel_for(ghost_pol, computeVelocity3D(var, rho, vel));
  Kokkos::fence();
  timers["calcSecond"].accumulate();

  for (int v = 0; v < cf.nv; ++v) {
    timers["flux"].reset();
    Kokkos::parallel_for(
        weno_pol, calculateWenoFluxesG(var, p, rho, vel, fluxx, fluxy, fluxz, cd, v));

    Kokkos::parallel_for(cell_pol,
                         advect3D(dvar, fluxx, fluxy, fluxz, v));
    Kokkos::fence();
    timers["flux"].accumulate();
  }

  timers["pressgrad"].reset();
  Kokkos::parallel_for(cell_pol, applyPressureGradient3D(dvar, p, cd));
  Kokkos::fence();
  timers["pressgrad"].accumulate();

  if (cf.buoyancy) {
    timers["buoyancy"].reset();
    Kokkos::parallel_for(cell_pol, computeBuoyancy3D(dvar, var, rho, cf.gAccel, cf.rhoRef));
    Kokkos::fence();
    timers["buoyancy"].accumulate();
  }

  if (cf.ceq) {
    timers["ceq"].reset();
    /**** CEQ ****/
    // find mac wavespeed
    Kokkos::parallel_reduce(cell_pol, maxWaveSpeed(var, p, rho, cd),
                            Kokkos::Max<double>(maxS));
#ifndef NOMPI
    MPI_Allreduce(&maxS, &maxS, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
#endif

    // calculate density and energy gradients and find ceq source terms
    Kokkos::parallel_for(cell_pol, calculateRhoGrad(var, rho, gradRho, cd));

    // update c equations
    Kokkos::parallel_for(cell_pol, updateCeq(dvar, var, gradRho, maxS, cd, cf.kap, cf.eps));

    /**** Apply CEQ ****/
     Kokkos::parallel_reduce(cell_pol, maxGradFunctor(var, cf.nv+0), Kokkos::Max<double>(maxC));
     Kokkos::parallel_reduce(cell_pol, maxGradFunctor(var, cf.nv+1), Kokkos::Max<double>(maxCh));
     Kokkos::parallel_reduce(cell_pol, maxGradFunctor(var, cf.nv+2), Kokkos::Max<double>(maxC1));
     Kokkos::parallel_reduce(cell_pol, maxGradFunctor(var, cf.nv+3), Kokkos::Max<double>(maxC2));
     Kokkos::parallel_reduce(cell_pol, maxGradFunctor(var, cf.nv+4), Kokkos::Max<double>(maxC3));
 #ifndef NOMPI
    MPI_Allreduce(&maxC,  &maxC,  1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxCh, &maxCh, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxC1, &maxC1, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxC2, &maxC2, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxC3, &maxC3, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
 #endif

    mu = maxC1;
    if (maxC2 > mu) mu = maxC2;
    if (maxC3 > mu) mu = maxC3;

    alpha = 0.0;
    beta = 0.0;
    betae = 0.0;

    double dxmag = sqrt(cf.dx*cf.dx + cf.dy*cf.dy + cf.dz*cf.dz);
    if (mu > 0 && maxCh > 0)
      alpha = (dxmag*dxmag / (mu * mu * maxCh)) * cf.alpha;
    if (maxC > 0) {
      beta = (dxmag*dxmag / maxC) * cf.beta;
      betae = (dxmag*dxmag / maxC) * cf.betae;
    }

    Kokkos::parallel_for(weno_pol, calculateCeqFlux(var, rho, mFlux, cFlux, cd));

    Kokkos::parallel_for(cell_pol, applyCeq(dvar, var, varx, vel, rho, mFlux, cFlux, cd, alpha, beta, betae));
    Kokkos::fence();
    timers["ceq"].accumulate();
  }
}
