#include "cart2d.hpp"
#include "vtk.hpp"

std::map<string,int> varxIds;

cart2d_func::cart2d_func(struct inputConfig &cf_, Kokkos::View<double *> &cd_)
    : rk_func(cf_, cd_) {

  int varxPtr = 0;

  grid  = FS4D("coords", cf.ni, cf.nj, cf.nk, 3);
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
    varxIds["Noise"]=varxPtr++;
    varxNames.push_back("Noise");
  }

  if (cf.particle == 1) {
    particles = FSP2D("particles", cf.p_np);         // Particle Array
    particlesH = Kokkos::create_mirror_view(particles);
  }

  assert(varNames.size() == cf.nvt);

  varxIds["X-Velocity"]=varxPtr++;
  varxNames.push_back("X-Velocity");

  varxIds["Y-Velocity"]=varxPtr++;
  varxNames.push_back("Y-Velocity");

  varxIds["Pressure"]=varxPtr++;
  varxNames.push_back("Pressure");

  varxIds["Temperature"]=varxPtr++;
  varxNames.push_back("Temperature");

  varxIds["Density"]=varxPtr++;
  varxNames.push_back("Density");

  for (int i=0; i<varxPtr; ++i)
    assert(varxIds[varxNames[i]]==i);
  for(std::pair<std::string,int> element : varxIds)
    assert(element.first.compare(varxNames[element.second]) == 0);
  assert(varxNames.size()==varxPtr);

  varx  = FS4D("varx", cf.ngi, cf.ngj, cf.ngk, varxPtr); // Extra Vars

  cd = mcd;


  timers["flux"] = fiestaTimer("Flux Calculation");
  timers["pressgrad"] = fiestaTimer("Pressure Gradient Calculation");
  timers["calcSecond"] = fiestaTimer("Secondary Variable Calculation");
  timers["solWrite"] = fiestaTimer("Solution Write Time");
  timers["resWrite"] = fiestaTimer("Restart Write Time");
  timers["statCheck"] = fiestaTimer("Status Check");
  timers["rk"] = fiestaTimer("Runge Stage Update");
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
  if (cf.particle == 1) {
    timers["padvect"] = fiestaTimer("Particle Advection");
    timers["pwrite"] = fiestaTimer("Particle Write");
    timers["psetup"] = fiestaTimer("Particle Setup");
  }
};

void cart2d_func::preStep() {}

void cart2d_func::postStep() {

      policy_f ghost_pol = policy_f({0, 0}, {cf.ngi, cf.ngj});
  if (cf.write_freq >0)
    if (cf.t % cf.write_freq == 0){
      //FS3D velx(varx,Kokkos::ALL,Kokkos::ALL,0,make_pair(0,cf.ndim));
      //FS2D px(varx,Kokkos::ALL,Kokkos::ALL,0,cf.ndim);
      //FS2D Tx(varx,Kokkos::ALL,Kokkos::ALL,0,cf.ndim+1);
      //FS2D rhox(varx,Kokkos::ALL,Kokkos::ALL,0,cf.ndim+2);
      

      // Calcualte Total Density and Pressure Fields
      timers["calcSecond"].reset();
      Kokkos::parallel_for(ghost_pol, calculateRhoPT2D(var, p, rho, T, cd));
      Kokkos::parallel_for(ghost_pol, computeVelocity2D(var, rho, vel));
      Kokkos::fence();
      timers["calcSecond"].accumulate();

      Kokkos::parallel_for( "Loop1", ghost_pol,
        KOKKOS_LAMBDA(const int i, const int j) {
          varx(i,j,0,varxIds["X-Velocity"]) = vel(i,j,0);
          varx(i,j,0,varxIds["Y-Velocity"]) = vel(i,j,1);
          varx(i,j,0,varxIds["Pressure"]) = p(i,j);
          varx(i,j,0,varxIds["Temperature"]) = T(i,j);
          varx(i,j,0,varxIds["Density"]) = rho(i,j);
      });
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
      double maxCh;
      Kokkos::parallel_reduce(cell_pol, maxCvar2D(var, 1, cd),
                              Kokkos::Max<double>(maxCh));
#ifndef NOMPI
      MPI_Allreduce(&maxCh, &maxCh, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
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

//#ifdef NOMPI
  if (cf.particle == 1) {
    // write particle data
#ifdef NOMPI
    if (cf.write_freq > 0) {
      if ((cf.t) % cf.write_freq == 0) {
        timers["pwrite"].reset();
        Kokkos::deep_copy(particlesH, particles);
        //        writeParticles(cf,particlesH);
        stringstream ss;
        ss << "particle-" << setw(7) << setfill('0') << cf.t << ".vtk";
        ofstream f;
        // f.open("particle.vtk");
        f.open(ss.str());
        f << "# vtk DataFile Version 4.2" << endl;
        f << "Test Particles" << endl;
        f << "ASCII" << endl;
        f << "DATASET POLYDATA" << endl;
        f << "POINTS " << cf.p_np << " float" << endl;
        for (int p = 0; p < cf.p_np; ++p) {
          f << particlesH(p).x << " " << particlesH(p).y << " "
            << "0.0" << endl;
        }
        f << "VERTICES " << cf.p_np << " " << cf.p_np * 2 << endl;
        for (int p = 0; p < cf.p_np; ++p) {
          f << "1 " << p << endl;
        }
        f << "POINT_DATA " << cf.p_np << endl;
        f << "SCALARS state float" << endl;
        f << "LOOKUP_TABLE default" << endl;
        for (int p = 0; p < cf.p_np; ++p) {
          f << particlesH(p).state << endl;
        }
        f << "SCALARS ci float" << endl;
        f << "LOOKUP_TABLE default" << endl;
        for (int p = 0; p < cf.p_np; ++p) {
          f << particlesH(p).ci << endl;
        }
        f << "SCALARS cj float" << endl;
        f << "LOOKUP_TABLE default" << endl;
        for (int p = 0; p < cf.p_np; ++p) {
          f << particlesH(p).cj << endl;
        }
        f.flush();
        f.close();

        Kokkos::fence();
        timers["pwrite"].accumulate();
      }
    } // end particle write
#endif

     // execution policy for all cells including ghost cells
     policy_f ghost_pol = policy_f({0, 0}, {cf.ngi, cf.ngj});
 
     // Calcualte Total Density and Pressure Fields
     timers["calcSecond"].reset();
     Kokkos::parallel_for(ghostPol, calculateRhoPT2D(var, p, rho, T, cd));
     Kokkos::fence();
     timers["calcSecond"].accumulate();
 
     // advect particles
     timers["padvect"].reset();
     Kokkos::parallel_for(
         cf.p_np, advectParticles2D(var, rho, grid, particles, cf.dt, cf.ng));
     Kokkos::fence();
     timers["padvect"].accumulate();
  }
//#endif
}
void cart2d_func::preSim() {

  if (cf.particle == 1) {
    timers["psetup"].reset();

    // find initial cell id
    policy_f grid_pol = policy_f({0, 0}, {cf.nci, cf.ncj});
    for (int p = 0; p < cf.p_np; ++p) {
      Kokkos::parallel_for(grid_pol,
                           findInitialCell2D(grid, particles, p, cf.ng));
    }
    Kokkos::fence();
    timers["psetup"].accumulate();

#ifdef NOMPI
    Kokkos::View<particleStruct2D *>::HostMirror particlesH =
        Kokkos::create_mirror_view(particles);
    Kokkos::deep_copy(particlesH, particles);

    // write particle data
    if (cf.write_freq > 0) {
      timers["pwrite"].reset();
      Kokkos::deep_copy(particlesH, particles);
      //   writeParticles(cf,particlesH);

      stringstream ss;
      ss << "particle-" << setw(7) << setfill('0') << cf.t << ".vtk";
      ofstream f;
      // f.open("particle.vtk");
      f.open(ss.str());
      f << "# vtk DataFile Version 4.2" << endl;
      f << "Test Particles" << endl;
      f << "ASCII" << endl;
      f << "DATASET POLYDATA" << endl;
      f << "POINTS " << cf.p_np << " float" << endl;
      for (int p = 0; p < cf.p_np; ++p) {
        f << particlesH(p).x << " " << particlesH(p).y << " "
          << "0.0" << endl;
      }
      f << "VERTICES " << cf.p_np << " " << cf.p_np * 2 << endl;
      for (int p = 0; p < cf.p_np; ++p) {
        f << "1 " << p << endl;
      }
      f << "POINT_DATA " << cf.p_np << endl;
      f << "SCALARS state float" << endl;
      f << "LOOKUP_TABLE default" << endl;
      for (int p = 0; p < cf.p_np; ++p) {
        f << particlesH(p).state << endl;
      }
      f << "SCALARS ci float" << endl;
      f << "LOOKUP_TABLE default" << endl;
      for (int p = 0; p < cf.p_np; ++p) {
        f << particlesH(p).ci << endl;
      }
      f << "SCALARS cj float" << endl;
      f << "LOOKUP_TABLE default" << endl;
      for (int p = 0; p < cf.p_np; ++p) {
        f << particlesH(p).cj << endl;
      }
      f.flush();
      f.close();

      Kokkos::fence();
    } // end initial write
#endif
  }
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
        face_pol, calculateStressTensor2dv(var, rho, T, stressx, stressy, cd));
    Kokkos::fence();
    timers["stress"].accumulate();

    timers["qflux"].reset();
    Kokkos::parallel_for(face_pol,
                         calculateHeatFlux2dv(var, rho, T, qx, qy, cd));
    Kokkos::fence();
    timers["qflux"].accumulate();

    timers["visc"].reset();
    Kokkos::parallel_for(cell_pol, applyViscousTerm2dv(dvar, var, rho, stressx,
                                                       stressy, qx, qy, cd));
    Kokkos::fence();
    timers["visc"].accumulate();
  }

  if (cf.ceq == 1) {
    double maxS;
    double maxG;
    double maxGh;
    double maxTau1;
    double maxTau2;
    double mu;
    double maxC;
    double maxCh;
    double maxT1;
    double maxT2;

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
    MPI_Allreduce(&maxS, &maxS, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxG, &maxG, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxGh, &maxGh, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxTau1, &maxTau1, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxTau1, &maxTau2, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxC, &maxC, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxCh, &maxCh, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxT1, &maxT1, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
    MPI_Allreduce(&maxT2, &maxT2, 1, MPI_DOUBLE, MPI_MAX, cf.comm);
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
