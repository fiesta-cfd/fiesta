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
#ifdef HAVE_MPI
#include "mpi.hpp"
#endif
#include <cassert>
#include "log2.hpp"

std::map<string,int> varxIds;

cart2d_func::cart2d_func(struct inputConfig &cf_) : rk_func(cf_) {
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
      varNames.push_back(cf.speciesName[v]+" Density");

  if (cf.visc) {
    qx      = FS2D("qx", cf.ngi, cf.ngj); // Heat Fluxes X direction
    qy      = FS2D("qy", cf.ngi, cf.ngj); // Heat Fluxes Y direction
    stressx = FS4D("stressx", cf.ngi, cf.ngj, 2, 2); // stress on x faces
    stressy = FS4D("stressy", cf.ngi, cf.ngj, 2, 2); // stress on y faces
  }

  if (cf.ceq) {
    gradRho = FS3D("gradRho", cf.ngi, cf.ngj, 4);    // Density Gradients
    m = FS5D("m",2,2,cf.ngi,cf.ngj,2);
    varNames.push_back("C");
    varNames.push_back("C_hat");
    varNames.push_back("Tau_1");
    varNames.push_back("Tau_2");
    varNames.push_back("Debug");
  }

  if (cf.noise) {
    noise = FS2D_I("noise", cf.ngi, cf.ngj);         // Noise Indicator
    //varxNames.push_back("Noise");
  }

  assert(varNames.size() == cf.nvt);

  varxNames.push_back("X-Velocity");
  varxNames.push_back("Y-Velocity");
  varxNames.push_back("Pressure");
  varxNames.push_back("Temperature");
  varxNames.push_back("Total Density");

  if (cf.ceq){
    varxNames.push_back("C_dvar_u");
    varxNames.push_back("C_dvar_v");
  }

  if (cf.noise){
    varxNames.push_back("Noise_I");
    varxNames.push_back("Noise_C");
    varxNames.push_back("Noise_D");
  }

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
  timers["advect"] = Timer::fiestaTimer("Advection Term Calculation");
  timers["pressgrad"] = Timer::fiestaTimer("Pressure Gradient Calculation");
  timers["calcSecond"] = Timer::fiestaTimer("Secondary Variable Calculation");
  timers["solWrite"] = Timer::fiestaTimer("Solution Write Time");
  timers["resWrite"] = Timer::fiestaTimer("Restart Write Time");
  timers["statCheck"] = Timer::fiestaTimer("Status Check");
  timers["rk"] = Timer::fiestaTimer("Runge Stage Update");
  timers["halo"] = Timer::fiestaTimer("Halo Exchanges");
  timers["bc"] = Timer::fiestaTimer("Boundary Conditions");
  if (cf.buoyancy) {
    timers["buoyancy"] = Timer::fiestaTimer("Buoyancy Term");
  }
  if (cf.visc) {
    timers["visc"] = Timer::fiestaTimer("Viscous Term Calculation");
  }
  if (cf.ceq) {
    timers["ceq"] = Timer::fiestaTimer("C-Equation");
  }
  if (cf.noise) {
    timers["noise"] = Timer::fiestaTimer("Noise Removal");
  }
};

void cart2d_func::preSim() {
  timers["calcSecond"].reset();
  policy_f ghost_pol = policy_f({0, 0}, {cf.ngi, cf.ngj});
  Kokkos::parallel_for(ghost_pol, calculateRhoPT2D(var, p, rho, T, cd));
  Kokkos::parallel_for(ghost_pol, computeVelocity2D(var, rho, vel));
  Kokkos::parallel_for(ghost_pol, copyExtraVars2D(varx, vel, p, rho, T));
  Kokkos::fence();
  timers["calcSecond"].accumulate();
}

void cart2d_func::preStep() {
}

void cart2d_func::compute() {
  policy_f ghost_pol = policy_f({0, 0}, {cf.ngi, cf.ngj});
  policy_f cell_pol  = policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});
  policy_f face_pol  = policy_f({cf.ng - 1, cf.ng - 1}, {cf.ngi - cf.ng, cf.ngj - cf.ng});

  // Calcualte Total Density and Pressure Fields
  timers["calcSecond"].reset();
  Kokkos::parallel_for(ghost_pol, calculateRhoPT2D(var, p, rho, T, cd));
  Kokkos::parallel_for(ghost_pol, computeVelocity2D(var, rho, vel));
  Kokkos::fence();
  timers["calcSecond"].accumulate();

  // Calculate and apply weno fluxes for each variable
  for (int v = 0; v < cf.nv; ++v) {
    timers["advect"].reset();
    if (cf.scheme == 3) {
      Kokkos::parallel_for( face_pol, computeFluxQuick2D(var, p, rho, fluxx, fluxy, cd, v));
    } else if (cf.scheme == 2) {
      Kokkos::parallel_for( face_pol, computeFluxCentered2D(var, p, rho, fluxx, fluxy, cd, v));
    } else {
      Kokkos::parallel_for( face_pol, computeFluxWeno2D(var, p, rho, vel, fluxx, fluxy, cf.dx, cf.dy, v));
    }
    //Kokkos::fence();
    Kokkos::parallel_for(cell_pol, advect2D(dvar, fluxx, fluxy, cd, v));
    Kokkos::fence();
    timers["advect"].accumulate();
  }

  // Apply Pressure Gradient Term
  timers["pressgrad"].reset();
  Kokkos::parallel_for(cell_pol, applyPressureGradient2D(dvar, p, cd));
  Kokkos::fence();
  timers["pressgrad"].accumulate();

  if (cf.buoyancy) {
    timers["buoyancy"].reset();
    Kokkos::parallel_for(cell_pol, computeBuoyancy(dvar, var, rho, cf.gAccel, cf.rhoRef));
    Kokkos::fence();
    timers["buoyancy"].accumulate();
  }

  if (cf.visc) {
    timers["visc"].reset();
    Kokkos::parallel_for(face_pol, calculateStressTensor2dv(var, rho, vel, stressx, stressy, cd));
    Kokkos::parallel_for(face_pol, calculateHeatFlux2dv(var, rho, T, qx, qy, cd));
    Kokkos::parallel_for(cell_pol, applyViscousTerm2dv(dvar, var, rho, vel, stressx, stressy, qx, qy, cd));
    Kokkos::fence();
    timers["visc"].accumulate();
  }

  if (cf.ceq) {
    FSCAL maxS,maxS_recv;
    FSCAL maxC,maxC_recv;
    FSCAL maxCh,maxCh_recv;

    timers["ceq"].reset();

    Kokkos::parallel_for(cell_pol, calculateRhoGrad2D(var, vel, rho, gradRho, cf.dx, cf.dy));

    Kokkos::parallel_reduce(cell_pol, maxWaveSpeed2D(var, p, rho, cd), Kokkos::Max<FSCAL>(maxS));
    Kokkos::parallel_reduce(cell_pol, maxCvar2D(var, 0, cd), Kokkos::Max<FSCAL>(maxC));
    Kokkos::parallel_reduce(cell_pol, maxCvar2D(var, 1, cd), Kokkos::Max<FSCAL>(maxCh));
    Kokkos::fence();

#ifdef HAVE_MPI
    MPI_Allreduce(&maxS, &maxS_recv, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxC, &maxC_recv, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxCh, &maxCh_recv, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
#endif

    maxS=maxS_recv;
    maxC=maxC_recv;
    maxCh=maxCh_recv;

    FSCAL alpha = (cf.dx*cf.dx + cf.dy*cf.dy)/(maxCh+1.0e-6)*cf.alpha;

    Kokkos::parallel_for(cell_pol, updateCeq2D(dvar, var, gradRho, maxS, cd, cf.kap, cf.eps));
    Kokkos::parallel_for(face_pol, computeCeqFlux2D(var, m, rho, alpha, cf.nv, maxCh));
    Kokkos::parallel_for(face_pol, computeCeqFaces2D(m, vel, cd));
    Kokkos::parallel_for(cell_pol, applyCeq2D(dvar, varx, m, cf.dx, cf.dy));

    Kokkos::fence();
    timers["ceq"].accumulate();
  }
}

void cart2d_func::postStep() {
  bool varsxNeeded=false;

  if ( (cf.write_freq >0) && ((cf.t+1) % cf.write_freq == 0) ) varsxNeeded=true;
  if ( (cf.stat_freq  >0) && ((cf.t+1) % cf.stat_freq  == 0) ) varsxNeeded=true;
  if (cf.ioThisStep) varsxNeeded = true;

  // compute secondary variables
  if (varsxNeeded){
    timers["calcSecond"].reset();
    policy_f ghost_pol = policy_f({0, 0}, {cf.ngi, cf.ngj});
    Kokkos::parallel_for(ghost_pol, calculateRhoPT2D(var, p, rho, T, cd));
    Kokkos::parallel_for(ghost_pol, computeVelocity2D(var, rho, vel));
    Kokkos::parallel_for(ghost_pol, copyExtraVars2D(varx, vel, p, rho, T));
    Kokkos::fence();
    timers["calcSecond"].accumulate();
  }

  if (cf.noise == 1) {
    int M = 0;
    int N = 0;
    FSCAL coff;

    if ((cf.nci + 1) % 2 == 0)
      M = (cf.nci + 1) / 2;
    else
      M = cf.nci / 2;

    if ((cf.ncj + 1) % 2 == 0)
      N = (cf.ncj + 1) / 2;
    else
      N = cf.ncj / 2;

    policy_f noise_pol = policy_f({0, 0}, {M+1, N+1});
    policy_f cell_pol = policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});

    if (cf.ceq == 1) {
      FSCAL maxCh,maxCh_recv;
      Kokkos::parallel_reduce(cell_pol, maxCvar2D(var, 1, cd), Kokkos::Max<FSCAL>(maxCh));
      #ifdef HAVE_MPI
        MPI_Allreduce(&maxCh, &maxCh_recv, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
      #endif
      maxCh=maxCh_recv;
      coff = cf.n_coff * maxCh;
    } else {
      coff = 0.0;
    }

    vector<int> noise_variables;
    if (cf.n_mode==1){ //velocity
      noise_variables.push_back(0);
      noise_variables.push_back(1);
    }else{
      noise_variables.push_back(2);
    }

    timers["noise"].reset();
    for (auto v : noise_variables) {
      Kokkos::parallel_for(noise_pol, detectNoise2D(var, varx, noise, cf.n_dh, coff, cd, v));
      for (int tau = 0; tau < cf.n_nt; ++tau) {
        Kokkos::parallel_for(cell_pol, removeNoise2D(dvar, var, varx, noise, cf.n_eta, cd, v));
      }
    }
    Kokkos::fence();
    timers["noise"].accumulate();
  }
}

void cart2d_func::postSim() {
}
